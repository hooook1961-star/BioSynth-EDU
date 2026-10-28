"""RDKit ligand depiction with atom/group-anchored residue interaction lines."""
from html import escape
import io
import math

from core.complex_analysis import ComplexAnalysisError, KINDS


def validated_ligand(result):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(result['smiles'])
    if mol is None or not mol.GetNumAtoms():
        raise ComplexAnalysisError('Не удалось восстановить химическую структуру лиганда для рисунка. Проверьте порядки связей в исходном комплексе.')
    heavy = {k: v for k, v in result['atoms'].items() if v['element'] not in {'H', 'D'}}
    mapping = {k: v for k, v in result['mapping'].items() if k in heavy}
    # Never connect a residue to an invented atom index. Validate full coverage
    # and element identities before using PLIP's canonical-SMILES atom map.
    if set(mapping) != set(heavy) or set(mapping.values()) != set(range(mol.GetNumAtoms())):
        raise ComplexAnalysisError('Нет однозначного соответствия атомов PDB и 2D-структуры. Таблица доступна, диаграмма не построена.')
    for serial, index in mapping.items():
        if mol.GetAtomWithIdx(index).GetSymbol().upper() != heavy[serial]['element']:
            raise ComplexAnalysisError('Нумерация атомов рисунка не совпала с PDB. Диаграмма не построена.')
    for contact in result['contacts']:
        if any(serial not in mapping for serial in contact['ligand_atoms']):
            raise ComplexAnalysisError('В отчёте есть взаимодействие с неопознанным атомом лиганда.')
    return mol


def diagram_layout(contacts):
    residues = sorted({(c['residue'], c['chain'], c['number']) for c in contacts}, key=lambda x: (x[1], x[2], x[0]))
    if len(residues) > 40:
        raise ComplexAnalysisError('На рисунке слишком много остатков. Отфильтруйте типы взаимодействий; полная таблица сохранена.')
    height = max(740, math.ceil(len(residues) / 2) * 64 + 240)
    positions = {}
    for i, residue in enumerate(residues):
        side, row = i % 2, i // 2
        positions[residue] = (1090 if side else 110, 150 + row * 64)
    return 1200, height, positions


def render_diagram(result, included_kinds=None, *, output_png=True):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.Mol(validated_ligand(result))
    rdDepictor.Compute2DCoords(mol)  # Copy only: deposited 3D pose is untouched.
    contacts = [c for c in result['contacts'] if included_kinds is None or c['kind'] in included_kinds]
    width, height, positions = diagram_layout(contacts)
    panel_width, panel_height = 700, height - 220
    dx, dy = 250, 90
    drawers = [rdMolDraw2D.MolDraw2DSVG(panel_width, panel_height)]
    if output_png:
        drawers.append(rdMolDraw2D.MolDraw2DCairo(panel_width, panel_height))
    for drawer in drawers:
        drawer.drawOptions().padding = 0.12
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol)
    anchors = {serial: drawers[0].GetDrawCoords(index) for serial, index in result['mapping'].items()
               if serial in result['atoms'] and result['atoms'][serial]['element'] not in {'H', 'D'}}
    for drawer in drawers:
        drawer.FinishDrawing()
    molecule_svg = drawers[0].GetDrawingText()
    molecule_svg = molecule_svg[molecule_svg.index('<svg'):]
    if output_png:
        from PIL import Image, ImageDraw, ImageFont
        image = Image.new('RGB', (width, height), 'white')
        molecule_png = Image.open(io.BytesIO(drawers[1].GetDrawingText())).convert('RGBA')
        image.paste(molecule_png, (dx, dy), molecule_png)
        pen = ImageDraw.Draw(image)
        try:
            font = ImageFont.truetype('DejaVuSans.ttf', 16)
            title_font = ImageFont.truetype('DejaVuSans.ttf', 23)
        except OSError:
            font = ImageFont.load_default(size=16)
            title_font = ImageFont.load_default(size=23)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<g transform="translate({dx},{dy})">{molecule_svg}</g>']

    def text(x, y, value, size=16, color='#172033'):
        svg.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans,Arial,sans-serif" font-size="{size}" fill="{color}">{escape(value)}</text>')
        if output_png:
            pen.text((x, y - size), value, fill=color, font=title_font if size > 20 else font)

    def line(x1, y1, x2, y2, color, dashed):
        dash = ' stroke-dasharray="6 5"' if dashed else ''
        svg.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="2"{dash}/>')
        if not output_png:
            return
        if dashed:
            length = math.hypot(x2 - x1, y2 - y1)
            for start in range(0, max(1, int(length)), 11):
                a, b = start / max(length, 1), min(start + 6, length) / max(length, 1)
                pen.line((x1+(x2-x1)*a, y1+(y2-y1)*a, x1+(x2-x1)*b, y1+(y2-y1)*b), fill=color, width=2)
        else:
            pen.line((x1, y1, x2, y2), fill=color, width=2)

    text(30, 40, 'Взаимодействия лиганда с белком', 24)
    text(30, 67, f"{result['ligand']} · PLIP {result['plip_version']} · расстояния в Å")
    for index, c in enumerate(contacts, 1):
        pts = [anchors[serial] for serial in c['ligand_atoms']]
        x = dx + sum(p.x for p in pts) / len(pts)
        y = dy + sum(p.y for p in pts) / len(pts)
        rx, ry = positions[(c['residue'], c['chain'], c['number'])]
        sx = rx + 90 if rx < width / 2 else rx - 90
        color = KINDS[c['kind']][3]
        line(sx, ry, x, y, color, c['kind'] != 'hydrogen_bond')
        text(sx + (12 if sx < width / 2 else -120), ry - 9 - (index % 3) * 15,
             f"{c['distance']:.2f} Å", color=color)
    for (name, chain, number), (x, y) in positions.items():
        svg.append(f'<rect x="{x-90}" y="{y-20}" width="180" height="40" rx="10" fill="#f1f5f9" stroke="#cbd5e1"/>')
        if output_png:
            pen.rounded_rectangle((x-90, y-20, x+90, y+20), radius=10, fill='#f1f5f9', outline='#cbd5e1')
        text(x-80, y+6, f'{name} {number} · {chain or "—"}')
    for i, (kind, (label, _, _, color)) in enumerate(KINDS.items()):
        x, y = 30 + (i % 3) * 390, height - 90 + (i // 3) * 30
        line(x, y, x+30, y, color, kind != 'hydrogen_bond')
        text(x+40, y+5, label)
    if not contacts:
        text(350, height-130, 'Для выбранных типов взаимодействий контактов нет.')
    text(30, height-12, 'Геометрические контакты выбранной позы; не оценка энергии связывания.')
    svg.append('</svg>')
    png_bytes = None
    if output_png:
        png = io.BytesIO()
        image.save(png, format='PNG')
        png_bytes = png.getvalue()
    return ''.join(svg).encode('utf-8'), png_bytes


def render_available_formats(result, included_kinds=None):
    """SVG is independent of Cairo, Pillow and installed PNG text fonts."""
    svg, _ = render_diagram(result, included_kinds, output_png=False)
    try:
        _, png = render_diagram(result, included_kinds, output_png=True)
        return svg, png, None
    except Exception as exc:
        import logging
        logging.getLogger(__name__).exception('PNG export failed; SVG is available')
        return svg, None, f'{type(exc).__name__}: {str(exc)[:1500]}'
