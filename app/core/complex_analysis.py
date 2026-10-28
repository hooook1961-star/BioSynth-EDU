"""Single-pose PDB validation and isolated, offline PLIP analysis.

No network requests, global PLIP configuration or shared upload files.
"""
import csv
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

MAX_BYTES = 10 * 1024 * 1024
AMINO_ACIDS = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE SEC PYL HIE HID HIP CYX'.split())
ADDITIVES = set('HOH WAT DOD GOL EDO PEG PGE PG4 PG5 PG6 MPD DMS BME EOH MOH IPA ACT ACE FMT TRS MES HEP EPE MOP BIC BIS CIT TAR TLA SO4 PO4 NO3 CO3 CL BR I NA K CA MG ZN MN FE CU CO NI CD HG NH4'.split())
KINDS = {
    'hydrophobic_interaction': ('Гидрофобный контакт', 'dist', 'Между атомами углерода', '#64748b'),
    'hydrogen_bond': ('Водородная связь', 'dist_d-a', 'Донор — акцептор', '#2563eb'),
    'salt_bridge': ('Солевой мостик', 'dist', 'Между центрами заряженных групп', '#dc2626'),
    'pi_stack': ('π–π взаимодействие', 'centdist', 'Между центрами колец', '#16a34a'),
    'pi_cation_interaction': ('Катион–π взаимодействие', 'dist', 'Заряженная группа — центр кольца', '#9333ea'),
}


class ComplexAnalysisError(ValueError):
    pass


def inspect_pdb(data):
    if not data or len(data) > MAX_BYTES:
        raise ComplexAnalysisError('Загрузите непустой PDB-файл размером не более 10 МБ.')
    try:
        text = data.decode('utf-8-sig')
    except UnicodeError as exc:
        raise ComplexAnalysisError('Не удалось прочитать текстовый PDB-файл.') from exc
    atoms, groups, models, protein = {}, {}, 0, False
    for line in text.splitlines():
        record = line[:6].strip()
        if record == 'MODEL':
            models += 1
        if record not in {'ATOM', 'HETATM'}:
            continue
        try:
            serial = int(line[6:11])
            coords = tuple(float(line[i:i+8]) for i in (30, 38, 46))
            resname, chain, number = line[17:20].strip(), line[21:22].strip(), int(line[22:26])
            element = line[76:78].strip().upper() or line[12:16].strip().lstrip('0123456789')[:1].upper()
            if not resname or serial <= 0 or not all(math.isfinite(v) for v in coords):
                raise ValueError()
        except (ValueError, IndexError) as exc:
            raise ComplexAnalysisError('Повреждена запись атома или его координаты в PDB.') from exc
        if serial in atoms:
            raise ComplexAnalysisError('Повторяются номера атомов. Загрузите комплекс с одной моделью и уникальной нумерацией.')
        if line[16:17].strip():
            raise ComplexAnalysisError('Обнаружены альтернативные положения атомов. Оставьте одну конформацию перед загрузкой.')
        if line[26:27].strip():
            raise ComplexAnalysisError('Коды вставки в номерах остатков пока не поддерживаются. Перенумеруйте остатки без кодов вставки.')
        key = f'{resname}:{chain}:{number}'
        atom = dict(serial=serial, element=element, coordinates=coords, residue=key)
        atoms[serial] = atom
        groups.setdefault(key, []).append(atom)
        protein = protein or (record == 'ATOM' and resname in AMINO_ACIDS)
    if models > 1:
        raise ComplexAnalysisError('Загрузите PDB-комплекс с одной моделью и одной позой лиганда.')
    if len(atoms) > 80000:
        raise ComplexAnalysisError('Комплекс слишком большой: допускается до 80 000 атомов.')
    if not protein:
        raise ComplexAnalysisError('В файле не найден белок. Нужен комплекс белка с лигандом.')
    candidates = []
    for key, members in groups.items():
        name = key.split(':')[0]
        heavy = [a for a in members if a['element'] not in {'H', 'D'}]
        if name not in AMINO_ACIDS | ADDITIVES and len(heavy) >= 2 and any(a['element'] == 'C' for a in heavy):
            if len(heavy) > 300:
                continue
            candidates.append(dict(key=key, atom_count=len(heavy)))
    if not candidates:
        raise ComplexAnalysisError('Не найден подходящий малый органический лиганд. Вода, ионы и распространённые добавки исключены.')
    return dict(text=text, atoms=atoms, candidates=candidates)


def parse_plip_report(xml_text, selected_key):
    """Keep atom IDs and distance definitions alongside display columns."""
    try:
        root = ET.fromstring(xml_text)
        site = next((s for s in root.findall('bindingsite') if ':'.join(
            (s.findtext('identifiers/' + tag) or '').strip() for tag in ('hetid', 'chain', 'position')) == selected_key), None)
        if site is None:
            raise ComplexAnalysisError('PLIP не распознал выбранный лиганд. Проверьте записи HETATM и химическую структуру комплекса.')
        if site.findtext('identifiers/composite') == 'True':
            raise ComplexAnalysisError('PLIP объединил несколько компонентов в составной лиганд. Для этой версии нужен один остаток лиганда.')
        mapping = {}
        raw_mapping = site.findtext('mappings/smiles_to_pdb') or ''
        for pair in raw_mapping.split(','):
            if pair.strip():
                smiles_index, pdb_index = map(int, pair.split(':'))
                mapping[pdb_index] = smiles_index - 1
        contacts = []
        for kind, (label, distance_field, definition, _) in KINDS.items():
            for node in site.findall(f'interactions/*/{kind}'):
                distance = float(node.findtext(distance_field))
                if not math.isfinite(distance) or distance <= 0:
                    raise ValueError('Invalid distance')
                if kind == 'hydrophobic_interaction':
                    ids = [int(node.findtext('ligcarbonidx'))]
                elif kind == 'hydrogen_bond':
                    field = 'acceptoridx' if node.findtext('protisdon') == 'True' else 'donoridx'
                    ids = [int(node.findtext(field))]
                else:
                    ids = [int(n.text) for n in node.findall('lig_idx_list/idx')]
                if not ids:
                    raise ValueError('Missing ligand atoms')
                contacts.append(dict(residue=node.findtext('restype') or '',
                                     chain=node.findtext('reschain') or '',
                                     number=int(node.findtext('resnr')), kind=kind,
                                     interaction=label, distance=distance,
                                     distance_definition=definition, ligand_atoms=ids))
        return dict(ligand=selected_key, smiles=site.findtext('identifiers/smiles') or '',
                    mapping=mapping, contacts=contacts,
                    plip_version=root.findtext('plipversion') or 'unknown')
    except ComplexAnalysisError:
        raise
    except (ET.ParseError, ValueError, TypeError) as exc:
        raise ComplexAnalysisError('Не удалось прочитать отчёт PLIP. Координаты и взаимодействия не подставлялись.') from exc


def restore_pdb_atom_ids(result, inspected, selected_key):
    """Translate PLIP 3's report IDs to input serials, without guessing offsets.

    PDBParser.proteinmap counts ATOM/HETATM rows from 1, with one skipped
    number after TER, even with NOFIX. Its report IDs are not necessarily the
    serial numbers in columns 7-11. Apply the same explicit map to contacts
    AND canonical SMILES mappings (including when IDs happen to overlap).
    See PLIP v3.0.1 plip/structure/preparation.py:PDBParser.parse_pdb.
    """
    report_to_input = {}
    report_id, after_ter = 0, False
    for line in inspected['text'].splitlines():
        record = line[:6].strip()
        if record == 'TER':
            after_ter = True
        elif record in {'ATOM', 'HETATM'}:
            report_id += 2 if after_ter else 1
            report_to_input[report_id] = int(line[6:11])
            after_ter = False
    ligand_atoms = {serial: atom for serial, atom in inspected['atoms'].items()
                    if atom['residue'] == selected_key}

    def original_id(index):
        serial = report_to_input.get(index)
        if serial not in ligand_atoms:
            raise ComplexAnalysisError('Не удалось однозначно сопоставить атомы отчёта PLIP с выбранным лигандом. Анализ остановлен.')
        return serial

    # Validate everything before returning a translated result. Never accept an
    # ID merely because it coincides with an input ligand serial number.
    contacts = [dict(c, ligand_atoms=[original_id(i) for i in c['ligand_atoms']])
                for c in result['contacts']]
    mapping = {original_id(i): index for i, index in result['mapping'].items()}
    return dict(result, contacts=contacts, mapping=mapping, atoms=ligand_atoms)


def analyze_complex(data, selected_key):
    inspected = inspect_pdb(data)
    if selected_key not in {c['key'] for c in inspected['candidates']}:
        raise ComplexAnalysisError('Выбранный лиганд отсутствует в загруженном комплексе.')
    with tempfile.TemporaryDirectory(prefix='biosynth-complex-') as folder:
        source = Path(folder) / 'complex.pdb'
        source.write_text(inspected['text'], encoding='utf-8', newline='')
        command = [sys.executable, '-m', 'plip.plipcmd', '-f', str(source),
                   '-o', folder, '-x', '--nofix', '--maxthreads', '1']
        try:
            completed = subprocess.run(command, capture_output=True, text=True,
                                       timeout=90, cwd=folder, check=False)
        except subprocess.TimeoutExpired as exc:
            raise ComplexAnalysisError('Анализ превысил 90 секунд. Используйте комплекс меньшего размера.') from exc
        except OSError as exc:
            raise ComplexAnalysisError('Не удалось запустить модуль анализа на сервере.') from exc
        if completed.returncode:
            if 'No module named' in completed.stderr or 'ImportError' in completed.stderr:
                raise ComplexAnalysisError('Анализ пока недоступен на сервере: необходимо установить PLIP и Open Babel.')
            raise ComplexAnalysisError('PLIP не смог обработать комплекс. Проверьте формат PDB, связность и подготовку лиганда.')
        reports = list(Path(folder).glob('*.xml'))
        if len(reports) != 1:
            raise ComplexAnalysisError('PLIP завершился без отчёта о взаимодействиях.')
        result = parse_plip_report(reports[0].read_text(encoding='utf-8'), selected_key)
    return restore_pdb_atom_ids(result, inspected, selected_key)


def table_rows(result):
    return [{'Остаток': c['residue'], 'Цепь': c['chain'] or '—', 'Номер': c['number'],
             'Тип взаимодействия': c['interaction'], 'Расстояние, Å': c['distance'],
             'Определение расстояния': c['distance_definition'],
             'Атомы лиганда (PDB)': ', '.join(map(str, c['ligand_atoms']))}
            for c in result['contacts']]


def contacts_csv(result):
    rows = table_rows(result)
    headers = ['Остаток', 'Цепь', 'Номер', 'Тип взаимодействия', 'Расстояние, Å',
               'Определение расстояния', 'Атомы лиганда (PDB)']
    output = io.StringIO(newline='')
    writer = csv.DictWriter(output, fieldnames=headers)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: ("'" + v if isinstance(v, str) and v.startswith(('=', '+', '-', '@')) else v)
                         for k, v in row.items()})
    return output.getvalue().encode('utf-8-sig')
