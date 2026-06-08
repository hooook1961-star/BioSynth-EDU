# translations.py
LANGUAGES = {
    "Русский": {
        "nav_title": "🧪BioSynth-EDU",
        "lang_label": "🌐 Выберите язык",
        "select_placeholder": "-- Выберите из списка --",
        # Sidebar
        "sidebar_select_mol": "🧪 Выбор молекулы",
        "sidebar_kaz_cat": "🇰🇿 База KZ",
        "sidebar_kaz_label": "Отечественные препараты и кейсы:",
        "sidebar_world_cat": "🌍 Стандартные примеры",
        "sidebar_world_label": "Примеры лекарственных веществ:",
        "sidebar_manual": "✍️ Ввести SMILES",
        "sidebar_manual_label": "Или вставьте SMILES ниже:",
        # Tabs
        "tab_3d": "🔬 Структура и свойства",
        "tab_admet": "📊 ADMET Анализ",
        "tab_docking": "🧬 Докинг",
        "tab_edu": "📚 Обучение",
        "tab_project": "🚀 Исследовательский проект",
        # Categories (World Examples)
        "cat_analgesic": "Анальгетик",
        "cat_stimulant": "Стимулятор",
        "cat_antipyretic": "Жаропонижающее",
        "cat_nsaid": "НПВС",
        "cat_antibiotic": "Антибиотик",
        "cat_alkaloid": "Алкалоид",
        "cat_neuro": "Нейромедиатор",
        # Tab 1: 3D & Info
        "btn_build_3d": "🏗️ Построить 3D",
        "btn_optimize": "✨ Оптимизировать (MMFF94)",
        "spinner_optimize": "Минимизация энергии...",
        "style_label": "Стиль отображения",
        "style_surface": "Поверхность",
        "style_sphere": "Сферы",
        "style_stick": "Палочки",
        "style_line": "Линии",
        "download_sdf": "💾 Скачать структуру (SDF)",
        "download_help": "SDF файл содержит 3D координаты атомов для работы в проф. софте (AutoDock, PyMOL и др.)",
        "info_select_mol": "Выберите молекулу слева и нажмите 'Построить 3D'",
        "header_ref": "📚 Справочник",
        "metric_mw": "М. вес",
        "metric_logp": "LogP",
        "metric_rot_bonds": "Вращающихся связей",
        "spinner_chembl": "Запрос к ChEMBL...",
        "status_label": "Статус:",
        "status_approved": "Одобрено",
        "mechanism_label": "🔬 Механизм действия",
        "no_chembl": "Биологическая активность в ChEMBL не найдена",
        "ext_links": "🔗 Внешние базы:",
        "btn_pubchem": "Профиль в PubChem",
        "btn_chembl": "Данные ChEMBL (IC50/Ki)",
        "btn_similarity": "Найти похожие",
        "warn_no_pubchem": "Данные в PubChem не найдены",
        "physchem.title": "🔮 Расширенный физико-химический анализ",
        "physchem.info.choose_or_enter_structure": "⚠️ Выберите или введите структуру молекулы.",
        "physchem.tab.conformers": "📊 Конформационный анализ",
        "physchem.tab.charges": "⚡ Распределение зарядов",
        "physchem.tab.descriptors": "📝 Топологические дескрипторы",
        "physchem.conformers.description": "**Анализ энергетического профиля конформеров (MMFF94; если параметры недоступны — UFF)**",
        "physchem.conformers.num_conformers": "Число конформеров",
        "physchem.conformers.calculate": "Рассчитать конформационный ансамбль",
        "physchem.conformers.spinner": "Поиск конформационных минимумов...",
        "physchem.conformers.show_energy_table": "Показать таблицу энергий",
        "physchem.conformers.download_sdf": "📥 Скачать пространственную структуру минимума (SDF)",
        "physchem.conformers.failed": "Не удалось построить или оптимизировать конформеры для этой структуры. Чаще всего это происходит для некорректного SMILES, солей, металлокомплексов или молекул без параметров силового поля.",
        "physchem.conformers.caption": "Нажмите кнопку, чтобы рассчитать конформационный ансамбль.",
        "physchem.charges.description": "**Картирование парциальных зарядов по методу Гастейгера-Марсили**",
        "physchem.charges.calculate": "Рассчитать и показать заряды",
        "physchem.charges.spinner": "Расчет 3D-геометрии и зарядов...",
        "physchem.charges.heavy_atoms_title": "**Заряды тяжёлых атомов**",
        "physchem.charges.failed": "Не удалось рассчитать заряды или построить 3D-структуру для этой молекулы.",
        "physchem.charges.caption": "Нажмите кнопку, чтобы рассчитать распределение зарядов.",
        "physchem.descriptors.description": "**Дескрипторы реакционной способности и полярной поверхности**",
        "physchem.descriptors.failed": "Не удалось рассчитать дескрипторы для этой структуры.",
        "physchem.table.conformer": "Конформер",
        "physchem.table.conformer_number": "№{number}",
        "physchem.table.energy_kcal_mol": "Энергия (ккал/моль)",
        "physchem.table.delta_energy_kcal_mol": "ΔE от минимума (ккал/моль)",
        "physchem.table.atom_index": "Индекс",
        "physchem.table.atom": "Атом",
        "physchem.table.charge": "Заряд",
        "physchem.table.descriptor": "Дескриптор",
        "physchem.table.description": "Описание",
        "physchem.table.value": "Значение",
        "physchem.file.global_minimum_sdf": "global_minimum_structure.sdf",
        "physchem.descriptor.MolWt.name": "MolWt",
        "physchem.descriptor.MolWt.description": "Молекулярная масса",
        "physchem.descriptor.ExactMolWt.name": "ExactMolWt",
        "physchem.descriptor.ExactMolWt.description": "Точная молекулярная масса",
        "physchem.descriptor.LogP.name": "LogP",
        "physchem.descriptor.LogP.description": "Crippen logP",
        "physchem.descriptor.MR.name": "MR",
        "physchem.descriptor.MR.description": "Молярная рефракция",
        "physchem.descriptor.TPSA.name": "TPSA",
        "physchem.descriptor.TPSA.description": "Топологическая полярная поверхность",
        "physchem.descriptor.HBA.name": "HBA",
        "physchem.descriptor.HBA.description": "Акцепторы H-связей",
        "physchem.descriptor.HBD.name": "HBD",
        "physchem.descriptor.HBD.description": "Доноры H-связей",
        "physchem.descriptor.RotatableBonds.name": "RotatableBonds",
        "physchem.descriptor.RotatableBonds.description": "Вращающиеся связи",
        "physchem.descriptor.RingCount.name": "RingCount",
        "physchem.descriptor.RingCount.description": "Число колец",
        "physchem.descriptor.AromaticRings.name": "AromaticRings",
        "physchem.descriptor.AromaticRings.description": "Ароматические кольца",
        "physchem.descriptor.HeavyAtoms.name": "HeavyAtoms",
        "physchem.descriptor.HeavyAtoms.description": "Тяжёлые атомы",
        "physchem.descriptor.FractionCSP3.name": "FractionCSP3",
        "physchem.descriptor.FractionCSP3.description": "Доля sp3-углерода",
        "physchem.descriptor.QED.name": "QED",
        "physchem.descriptor.QED.description": "Лекарственное подобие QED",
        "admet_header": "📊 ADMET Анализ и интерпретация",
        "admet_instructions": """
        Для проведения анализа ADMET:
        1. Нажмите кнопку **«Открыть ADMETlab 3.0»** ниже.
        2. Скопируйте ваш SMILES и вставьте его в поле ввода на сайте.
        3. После завершения расчета скачайте результат в формате **CSV**.
        4. Загрузите файл с помощью кнопки Upload для автоматической интерпретации.
        """,
        "btn_open_admetlab": "🌐 Открыть ADMETlab 3.0",
        "btn_open_swissadme": "🧪 Альтернатива: SwissADME",
        "uploader_label": "📥 Загрузите CSV файл с результатами ADMETlab",
        "expander_raw": "👁️ Посмотреть все сырые данные CSV",
        "subheader_interp": "📝 Краткая интерпретация ключевых параметров",
        "metric_lipophilicity": "LogP (Липофильность)",
        "metric_bbb": "BBB (Проницаемость)",
        "metric_tox": "Токсичность / hERG",
        "status_optimal": "✅ Оптимально.",
        "status_extreme": "⚠️ Крайние значения.",
        "status_bbb_yes": "🧠 Проникает через ГЭБ.",
        "status_bbb_no": "🛡️ Низкий риск для ЦНС.",
        "status_tox_high": "💔 Высокий риск.",
        "status_tox_low": "❤️ Риск низкий.",
        "header_lipinski": "🧐 Соответствие правилам Drug-like",
        "lipinski_success": "🌟 Молекула полностью соответствует правилу Липинского!",
        "lipinski_warn": "⚠️ Нарушений правила Липинского: ",
        "lipinski_info": "💡 Напоминание: Допускается 1 нарушение для сохранения 'drug-likeness'.",
        "error_interp": "❌ Ошибка интерпретации: ",
                # Advanced CNS / BBB / QSAR module
        "advanced_cns_title": "🧠 Продвинутое исследование CNS / BBB / P-gp",
        "advanced_cns_caption": "Отдельный модуль для углублённого QSAR/ADMET-анализа молекул",
        "advanced_cns_body": (
            "Этот модуль предназначен для более подробного исследования молекулы по SMILES. "
            "Он позволяет не только получить ADMET-прогноз, но и разобрать, какие свойства "
            "молекулы влияют на проникновение через гематоэнцефалический барьер, риск P-gp "
            "эффлюкса и потенциальную CNS-доступность."
        ),
        "advanced_cns_features_title": "В модуле доступны:",
        "advanced_cns_features": [
            "расчёт RDKit-дескрипторов и BBB/ГЭБ-профиля;",
            "разбор влияния LogP, TPSA, MW, HBD/HBA, pKa, заряда и P-gp;",
            "матрица BBB × P-gp для объяснения пассивной проницаемости и активного эффлюкса;",
            "What-if лаборатория для учебного изменения дескрипторов;",
            "учебный отчёт по молекуле;",
            "дополнительный ML/SHAP-разбор RandomForest-моделей."
        ],
        "advanced_cns_note": (
            "Модуль предназначен для учебного и исследовательского QSAR/ADMET-анализа. "
            "Результаты являются in silico-прогнозом и не заменяют экспериментальную проверку."
        ),
        "advanced_cns_button": "Открыть модуль продвинутого исследования CNS",
        "docking_header": "🛠️ Подготовка лиганда к докингу",
        "docking_mol_ready": "✅ 3D-структура обнаружена и готова к обработке.",
        "docking_checklist_title": "**Чек-лист подготовки:**",
        "docking_checklist_items": """
        1. Добавление неявных водородов.
        2. Генерация 3D-конформации.
        3. Минимизация энергии (силовое поле MMFF94).
        4. **Определение активных торсионных углов (PDBQT).**
        """,
        "btn_run_prep": "⚙️ Запустить полную подготовку",
        "spinner_meeko": "Работают Meeko и RDKit: расчет зарядов и торсионов...",
        "docking_success_info": "Лиганд готов! Рассчитаны торсионы и заряды.",
        "docking_error": "Ошибка при подготовке PDBQT.",
        "docking_note_student": "ℹ️ **Заметка для студентов:** Докинг имитирует 'ключ и замок'. Чтобы ключ подошел, он должен иметь правильные углы связей.",
        "btn_download_pdbqt": "📥 Скачать готовый PDBQT",
        "docking_view_label": "🔍 **Просмотр подготовленной структуры:**",
        "docking_caption": "Оптимизированная 3D-модель (подготовлено для анализа).",
        "docking_next_steps_header": "🎓 Что дальше?",
        "docking_next_steps_text": """
        После подготовки лиганда вам необходимо:
        1. Подготовить **белок-мишень** (удалить воду, добавить заряды).
        2. Определить **Grid Box** (координаты активного центра).
        3. Запустить расчет в AutoDock Vina или аналогичном ПО.
        """,
        "docking_stage2_title": "🎯 Поиск мишени в базе scPDB",
        "docking_stage2_desc": "На основе структуры вашего лиганда система проведет мгновенный векторный QSAR-скрининг по базе данных scPDB.",
        "btn_run_screening": "🚀 Запустить ин silico поиск мишеней",
        "spinner_screening": "Проводится сравнение со структурами 17 188 лигандов scPDB...",
        "screening_success": "⚡ Скрининг успешно завершен!",
        "top_match_title": "🏆 **Наилучшее совпадение по сходству структуры:** Папка / ID сайта: **{pdb_id}**",
        "top_match_reason": "🧬 *Обоснование:* Индекс сходства Танимото с нативным лигандом: {sim:.2f}.",
        "top_match_instruction": "💡 *Инструкция для студента:* Ваша молекула наиболее похожа на нативный лиганд этого белка.",
        "table_title": "📊 **ТОП-15 потенциальных биологических мишеней для докинга:**",
        "col_pdb_id": "PDB ID / Сайт",
        "col_tanimoto": "Индекс Танимото (Сходство)",
        "col_pkd": "Прогноз pKd",
        "pdb_lookup_title": "🔍 Быстрый переход к структуре белка в RCSB PDB",
        "pdb_lookup_label": "Введите 4-значный PDB ID (например, 4E5H):",
        "btn_go_to_pdb": "Открыть карточку мишени {pdb_id} на сайте rcsb.org",
        "pdb_error_length": "⚠️ PDB ID должен состоять ровно из 4 символов.",
        "docking_warn_no_3d": "⚠️ Сначала постройте 3D модель на первой вкладке!",
            # scPDB target screening / ligand-based target hypotheses
        "docking_main_title": "🧬 Молекулярный докинг: Подбор сайта связывания",
        "target_screening_header": "🎯 Гипотезы мишеней по scPDB",
        "target_screening_subtitle": "Поиск возможных мишеней по сходству с нативными лигандами scPDB",
        "target_screening_input_label": "Введите SMILES для поиска target-гипотез",
        "target_screening_button": "Запустить поиск по scPDB",
        "target_screening_top15": "Топ-15 target-гипотез",
        "target_screening_best_match": "Лучшее совпадение",
        "target_screening_method": "Метод",
        "target_method_short": "scPDB ligand similarity screening",
        "target_method_note": (
            "Метод предлагает target-гипотезы на основе сходства введённой молекулы "
            "с нативными лигандами scPDB. Если доступен расширенный индекс, ranking "
            "дополнительно учитывает совместимость с признаками сайта связывания, cavity "
            "и interaction fingerprint. Результат является вычислительной гипотезой, "
            "а не docking-расчётом, прогнозом affinity или экспериментальным доказательством активности."
        ),
        "target_error_invalid_smiles": "Некорректный SMILES. Проверьте запись молекулы.",
        "target_error_db_unavailable": "База данных scPDB недоступна или пуста.",
        "target_error_unknown": "Ошибка при выполнении scPDB-скрининга.",
        "target_no_hits_message": (
            "Не найдено достаточно близких совпадений с нативными лигандами scPDB. "
            "Это не означает отсутствие мишеней, но ligand-based поиск не дал сильной гипотезы."
        ),
        "target_candidate_name": "Гипотеза мишени из scPDB: PDB {pdb_id}",
        "target_reason_ligand_similarity": (
            "Введённая молекула похожа на нативный лиганд комплекса {pdb_id}. "
            "Индекс Танимото = {similarity:.2f}; уровень: {similarity_label}."
        ),
        "target_student_interpretation": (
            "Это ligand-based гипотеза: если молекула похожа на лиганд, найденный "
            "в комплексе с белком, соответствующий белок можно рассматривать как "
            "возможную мишень для дальнейшего анализа."
        ),
        "target_limitation": (
            "Этот результат не является docking-расчётом, прогнозом affinity или "
            "экспериментальным доказательством активности."
        ),
        "target_similarity_high": "высокое сходство",
        "target_similarity_moderate": "умеренное сходство",
        "target_similarity_weak": "слабое сходство",
        "target_similarity_low": "низкое сходство",
        "target_col_rank": "Ранг",
        "target_col_pdb": "PDB ID",
        "target_col_similarity": "Tanimoto",
        "target_col_similarity_level": "Уровень сходства",
        "target_col_score": "Score, %",
        "target_metric_best_similarity": "Лучшее Tanimoto-сходство",
        "target_metric_hits": "Совпадений выше порога",
        "target_metric_database_size": "Размер базы scPDB",
        "target_min_similarity": "Минимальный порог similarity",
        "target_best_match_header": "Лучшее совпадение",
        "target_screening_stats_header": "Сводка scPDB-скрининга",
	"docking_stage3_title": "Docking box и конфигурация AutoDock Vina",
	"docking_stage3_desc": "Выберите одну из 15 scPDB-гипотез. Приложение покажет координаты активного сайта и подготовит конфигурационные файлы для AutoDock Vina.",
	"docking_select_candidate_label": "Выберите PDB / scPDB-запись для докинга",
	"docking_col_entry_id": "scPDB entry",
	"docking_col_box_center": "Центр бокса, Å",
	"docking_col_box_size": "Размер бокса, Å",
	"docking_col_box_status": "Статус бокса",
	"docking_metric_selected_pdb": "Выбранный PDB",
	"docking_metric_selected_entry": "scPDB entry",
	"docking_metric_box_status": "Docking box",
	"docking_box_coordinates_title": "#### Координаты docking box",
	"docking_coord_param": "Параметр",
	"docking_coord_value": "Значение",
	"docking_box_source_caption": "Центр рассчитан из: {center_source}. Размер рассчитан из: {size_source}.",
	"docking_box_unavailable": "Координаты docking box недоступны для этой записи.",
	"docking_config_title": "#### Готовый docking_config для AutoDock Vina",
	"docking_config_unavailable": "Конфигурационный файл пока недоступен для этой записи.",
	"btn_download_config_txt": "Скачать docking_config.txt",
	"btn_download_config_yml": "Скачать docking_config.yml",
	"docking_receptor_title": "#### Receptor / белковая структура",
	"docking_receptor_available": "Для этой записи доступен подготовленный receptor.pdbqt.",
	"docking_receptor_not_available": "Подготовленный receptor.pdbqt не хранится в облачной версии приложения.",
	"docking_receptor_instruction": "Скачайте структуру белка из RCSB PDB и подготовьте receptor.pdbqt локально в AutoDockTools, Open Babel или Meeko. В облаке BioSynth-EDU хранит только компактный scPDB-index, координаты docking box и конфигурации.",
	"docking_local_command_title": "#### Пример локального запуска AutoDock Vina",
	"docking_requires_molecule": "Сначала введите SMILES и подготовьте молекулу.",
        "project_warning": "👈 Пожалуйста, выберите соединение в боковой панели (База KZ), чтобы начать проект.",
        "project_mol_header": "Объект исследования",
        "project_smiles_label": "SMILES код",
        "project_task_header": "Задание для выполнения",
        "project_task_1": "1. Проведите прогноз биологической активности (PASS Online).",
        "project_task_2": "2. Оцените соответствие правилу Липинского (вкладка ADMET).",
        "project_task_3": "3. Подготовьте данные для итогового отчета.",
        "project_ketcher_header": "Редактор молекул (Ketcher)",
        "project_mod_warning": "Внимание: структура была изменена!",
        "project_btn_analyze": "🧪 Анализировать новую структуру",
        "tab_project": "🚀 Исследовательский проект",
        "project_desc": "**Исследовательский проект по органической химии для студентов и учащихся**",
        "project_start_hint": "👈 Для начала работы выберите молекулу в боковой панели\n(раздел «База kz» или поле ручного ввода SMILES)",
        "selected_mol": "Выбранное соединение",
        "input_struct": "🔬 Введённая структура",
        "current_smiles": "Текущий SMILES",
        "kz_info": "🇰🇿 Сведения о казахстанской разработке",
        "authors": "Авторы",
        "description": "Описание",
        "classification": "Классификация",
        "year": "Год",
        "info_missing": "Информация отсутствует",
        "task_title": "📝 Задание на исследовательский проект",
        "task_full": "Открыть полное задание",
        "task_step_1": "1. Проведите прогноз спектра биологической активности с помощью **PASS Online**.",
        "task_step_2": "2. Оцените фармакологические свойства (во вкладке **ADMET**) и проверьте выполнение правила Липинского.",
        "task_step_3": "3. В редакторе ниже **измените структуру молекулы** (добавьте/уберите функциональные группы, измените заместители) и проанализируйте, как это повлияло на свойства, повторив шаги 1 и 2.",
        "task_step_4": "4. Сформулируйте выводы и рекомендации для доклада.",
        "editor_title": "🧪 Редактор структуры молекулы",
        "editor_task": "Задание: Измените структуру ниже и нажмите «Применить изменения».",
        "structure_changed": "✅ Структура изменена в редакторе!",
        "new_smiles": "Новый SMILES",
        "apply_btn": "🔄 Применить изменения и обновить проект",
        "change_applied": "✅ Изменения применены!",
        "download_btn": "💾 Скачать изменённый SMILES",
        "change_editor_hint": "👆 Измените молекулу в редакторе выше",
        "lectures_header": "Видео-лекции по теме",
        "pass_lecture": "🎥 Лекция: PASS",
        "adme_lecture": "🎥 Лекция: ADME",
        "pubchem_lecture": "🎥 Лекция: PubChem",
        "tools_header": "Инструменты анализа",
        "pass_online": "🌐 PASS Online",
        "swissadme": "🧪 SwissADME",
        "pubchem_search": "📊 PubChem Search",
        "quiz_btn": "📝 Пройти тест и получить вопросы к защите",
        "quiz_title": "Интеллектуальный тренажер BioSynth-EDU",
        "quiz_part_1": "Часть 1: Тестирование",
        "quiz_part_2": "Часть 2: Вопросы для подготовки к докладу",
        "check_result": "Проверить результат",
        "quiz_score": "Ваш результат",
        "question_label": "Вопрос",
        "defense_questions": "Вопросы для подготовки к докладу",
        "close": "Закрыть",
        "mol_not_selected": "⚠️ Молекула не выбрана",
        "how_start": "Как начать:",
        "start_step_1": "• Выберите соединение из Казахстанского каталога в боковой панели",
        "start_step_2": "• Или введите SMILES вручную в боковой панели",
        "download_filename": "modified_molecule"
    },
    "Қазақша": {
        "nav_title": "BioSynth-EDU",
        "lang_label": "🌐 Тілді таңдаңыз",
        "select_placeholder": "-- Тізімнен таңдаңыз --",
        # Sidebar
        "sidebar_select_mol": "🧪 Молекуланы таңдау",
        "sidebar_kaz_cat": "🇰🇿 Қазақстандық әзірлемелер",
        "sidebar_kaz_label": "Отандық препараттар мен кейстер:",
        "sidebar_world_cat": "🌍 Стандартты мысалдар",
        "sidebar_world_label": "Дәрілік заттардың мысалдары:",
        "sidebar_manual": "✍️ SMILES енгізу",
        "sidebar_manual_label": "Немесе төмендегі SMILES-ті қойыңыз:",
        # Tabs
        "tab_3d": "🔬 Құрылымы мен қасиеттері",
        "tab_admet": "📊 ADMET Талдау",
        "tab_docking": "🧬 Докинг",
        "tab_edu": "📚 Оқыту",
        "tab_project": "🚀Зерттеу жобасы",
        # Caterogies (World Examples)
        "cat_analgesic": "Анальгетик",
        "cat_stimulant": "Стимулятор",
        "cat_antipyretic": "Ыстықты түсіретін",
        "cat_nsaid": "ҚАҚД",
        "cat_antibiotic": "Антибиотик",
        "cat_alkaloid": "Алкалоид",
        "cat_neuro": "Нейромедиатор",
        # Tab 1: 3D & Info
        "btn_build_3d": "🏗️ 3D құрастыру",
        "btn_optimize": "✨ Оңтайландыру (MMFF94)",
        "spinner_optimize": "Энергияны азайту...",
        "style_label": "Көрсету стилі",
        "style_surface": "Беттік",
        "style_sphere": "Сфералар",
        "style_stick": "Таяқшалар",
        "style_line": "Сызықтар",
        "download_sdf": "💾 Құрылымды жүктеу (SDF)",
        "download_help": "SDF файлында кәсіби бағдарламаларға (AutoDock, PyMOL т.б.) арналған атомдардың 3D координаттары бар",
        "info_select_mol": "Сол жақтан молекуланы таңдап, '3D құрастыру' батырмасын басыңыз",
        "header_ref": "📚 Анықтамалық",
        "metric_mw": "М. салмағы",
        "metric_logp": "LogP",
        "metric_rot_bonds": "Айналмалы байланыстар",
        "spinner_chembl": "ChEMBL-ге сұраныс жіберу...",
        "status_label": "Мәртебесі:",
        "status_approved": "Мақұлданған",
        "mechanism_label": "🔬 Әсер ету механизмі",
        "no_chembl": "ChEMBL-де биологиялық белсенділік табылмады",
        "ext_links": "🔗 Сыртқы дерекқорлар:",
        "btn_pubchem": "PubChem профилі",
        "btn_chembl": "ChEMBL деректері (IC50/Ki)",
        "btn_similarity": "Ұқсастарды табу",
        "warn_no_pubchem": "PubChem-нен деректер табылмады",
        "physchem.title": "🔮 Кеңейтілген физика-химиялық талдау",
        "physchem.info.choose_or_enter_structure": "⚠️ Молекула құрылымын таңдаңыз немесе енгізіңіз.",
        "physchem.tab.conformers": "📊 Конформациялық талдау",
        "physchem.tab.charges": "⚡ Зарядтардың таралуы",
        "physchem.tab.descriptors": "📝 Топологиялық дескрипторлар",
        "physchem.conformers.description": "**Конформерлердің энергетикалық профилін талдау (MMFF94; параметрлер қолжетімсіз болса — UFF)**",
        "physchem.conformers.num_conformers": "Конформерлер саны",
        "physchem.conformers.calculate": "Конформациялық ансамбльді есептеу",
        "physchem.conformers.spinner": "Конформациялық минимумдарды іздеу...",
        "physchem.conformers.show_energy_table": "Энергиялар кестесін көрсету",
        "physchem.conformers.download_sdf": "📥 Минимумның кеңістіктік құрылымын жүктеп алу (SDF)",
        "physchem.conformers.failed": "Бұл құрылым үшін конформерлерді құру немесе оңтайландыру мүмкін болмады. Бұл көбіне қате SMILES, тұздар, металл кешендері немесе күш өрісі параметрлері жоқ молекулалар үшін болады.",
        "physchem.conformers.caption": "Конформациялық ансамбльді есептеу үшін батырманы басыңыз.",
        "physchem.charges.description": "**Гастейгер-Марсили әдісі бойынша парциалдық зарядтарды картаға түсіру**",
        "physchem.charges.calculate": "Зарядтарды есептеу және көрсету",
        "physchem.charges.spinner": "3D-геометрия мен зарядтарды есептеу...",
        "physchem.charges.heavy_atoms_title": "**Ауыр атомдардың зарядтары**",
        "physchem.charges.failed": "Бұл молекула үшін зарядтарды есептеу немесе 3D-құрылым құру мүмкін болмады.",
        "physchem.charges.caption": "Зарядтардың таралуын есептеу үшін батырманы басыңыз.",
        "physchem.descriptors.description": "**Реакциялық қабілеттілік және полярлық бет дескрипторлары**",
        "physchem.descriptors.failed": "Бұл құрылым үшін дескрипторларды есептеу мүмкін болмады.",
        "physchem.table.conformer": "Конформер",
        "physchem.table.conformer_number": "№{number}",
        "physchem.table.energy_kcal_mol": "Энергия (ккал/моль)",
        "physchem.table.delta_energy_kcal_mol": "Минимумнан ΔE (ккал/моль)",
        "physchem.table.atom_index": "Индекс",
        "physchem.table.atom": "Атом",
        "physchem.table.charge": "Заряд",
        "physchem.table.descriptor": "Дескриптор",
        "physchem.table.description": "Сипаттама",
        "physchem.table.value": "Мәні",
        "physchem.file.global_minimum_sdf": "global_minimum_structure.sdf",
        "physchem.descriptor.MolWt.name": "MolWt",
        "physchem.descriptor.MolWt.description": "Молекулалық масса",
        "physchem.descriptor.ExactMolWt.name": "ExactMolWt",
        "physchem.descriptor.ExactMolWt.description": "Нақты молекулалық масса",
        "physchem.descriptor.LogP.name": "LogP",
        "physchem.descriptor.LogP.description": "Crippen logP",
        "physchem.descriptor.MR.name": "MR",
        "physchem.descriptor.MR.description": "Молярлық рефракция",
        "physchem.descriptor.TPSA.name": "TPSA",
        "physchem.descriptor.TPSA.description": "Топологиялық полярлық бет ауданы",
        "physchem.descriptor.HBA.name": "HBA",
        "physchem.descriptor.HBA.description": "H-байланыс акцепторлары",
        "physchem.descriptor.HBD.name": "HBD",
        "physchem.descriptor.HBD.description": "H-байланыс донорлары",
        "physchem.descriptor.RotatableBonds.name": "RotatableBonds",
        "physchem.descriptor.RotatableBonds.description": "Айналмалы байланыстар",
        "physchem.descriptor.RingCount.name": "RingCount",
        "physchem.descriptor.RingCount.description": "Сақиналар саны",
        "physchem.descriptor.AromaticRings.name": "AromaticRings",
        "physchem.descriptor.AromaticRings.description": "Ароматтық сақиналар",
        "physchem.descriptor.HeavyAtoms.name": "HeavyAtoms",
        "physchem.descriptor.HeavyAtoms.description": "Ауыр атомдар",
        "physchem.descriptor.FractionCSP3.name": "FractionCSP3",
        "physchem.descriptor.FractionCSP3.description": "sp3-көміртектің үлесі",
        "physchem.descriptor.QED.name": "QED",
        "physchem.descriptor.QED.description": "QED дәрілік ұқсастығы",
        "admet_header": "📊 ADMET талдау және интерпретация",
        "admet_instructions": """
        ADMET талдауын жүргізу үшін:
        1. Төмендегі **«ADMETlab 3.0 ашу»** батырмасын басыңыз.
        2. SMILES кодын көшіріп, сайттағы енгізу өрісіне қойыңыз.
        3. Есептеу аяқталғаннан кейін нәтижені **CSV** форматында жүктеңіз.
        4. Автоматты интерпретация үшін файлды Upload батырмасы арқылы жүктеңіз.
        """,
        "btn_open_admetlab": "🌐 ADMETlab 3.0 ашу",
        "btn_open_swissadme": "🧪 Балама: SwissADME",
        "uploader_label": "📥 ADMETlab нәтижелері бар CSV файлын жүктеңіз",
        "expander_raw": "👁️ Барлық шикі CSV деректерін көру",
        "subheader_interp": "📝 Негізгі параметрлердің қысқаша интерпретациясы",
        "metric_lipophilicity": "LogP (Липофильділік)",
        "metric_bbb": "BBB (Өтімділік)",
        "metric_tox": "Токсикология / hERG",
        "status_optimal": "✅ Оңтайлы.",
        "status_extreme": "⚠️ Шеткі мәндер.",
        "status_bbb_yes": "🧠 ГЭБ арқылы өтеді.",
        "status_bbb_no": "🛡️ ОЖЖ үшін қауіп төмен.",
        "status_tox_high": "💔 Жоғары қауіп.",
        "status_tox_low": "❤️ Қауіп төмен.",
        "header_lipinski": "🧐 Drug-like ережелеріне сәйкестігі",
        "lipinski_success": "🌟 Молекула Липинский ережесіне толық сәйкес келеді!",
        "lipinski_warn": "⚠️ Липинский ережесін бұзу саны: ",
        "lipinski_info": "💡 Ескерту: 'drug-likeness' қасиетін сақтау үшін 1 бұзушылыққа жол беріледі.",
        "error_interp": "❌ Интерпретация қатесі: ",
                # Advanced CNS / BBB / QSAR module
        "advanced_cns_title": "🧠 CNS / BBB / P-gp кеңейтілген зерттеу модулі",
        "advanced_cns_caption": "Молекулаларды терең QSAR/ADMET талдауға арналған жеке модуль",
        "advanced_cns_body": (
            "Бұл модуль SMILES бойынша молекуланы тереңірек зерттеуге арналған. "
            "Ол тек ADMET-болжам беріп қана қоймай, молекуланың қандай қасиеттері "
            "гематоэнцефалдық бөгеттен өтуге, P-gp эффлюкс қаупіне және ықтимал "
            "CNS-қолжетімділікке әсер ететінін түсіндіреді."
        ),
        "advanced_cns_features_title": "Модульде қолжетімді:",
        "advanced_cns_features": [
            "RDKit дескрипторларын және BBB профилін есептеу;",
            "LogP, TPSA, MW, HBD/HBA, pKa, заряд және P-gp әсерін талдау;",
            "пассивті өткізгіштік пен белсенді эффлюксті түсіндіретін BBB × P-gp матрицасы;",
            "дескрипторларды оқу мақсатында өзгертуге арналған What-if зертханасы;",
            "молекула бойынша оқу есебін қалыптастыру;",
            "RandomForest модельдері үшін қосымша ML/SHAP талдауы."
        ],
        "advanced_cns_note": (
            "Модуль оқу және зерттеу мақсатындағы QSAR/ADMET талдауға арналған. "
            "Нәтижелер in silico-болжам болып табылады және эксперименттік тексеруді алмастырмайды."
        ),
        "advanced_cns_button": "CNS кеңейтілген зерттеу модулін ашу",
        "docking_header": "🛠️ Лигандты докингке дайындау",
        "docking_mol_ready": "✅ 3D-құрылым табылды және өңдеуге дайын.",
        "docking_checklist_title": "**Дайындық тізімі (чек-лист):**",
        "docking_checklist_items": """
        1. Жасырын сутегілерді қосу.
        2. 3D-конформацияны генерациялау.
        3. Энергияны минимизациялау (MMFF94 күш өрісі).
        4. **Белсенді торсиондық бұрыштарды анықтау (PDBQT).**
        """,
        "docking_stage2_title": "🎯 scPDB дерекқорынан нысананы іздеу",
        "docking_stage2_desc": "Лигандтың құрылымы негізінде жүйе scPDB дерекқоры бойынша лезде векторлық QSAR-скрининг жүргізеді.",
        "btn_run_screening": "🚀 Нысаналарды ин silico іздеуді іске қосу",
        "spinner_screening": "scPDB дерекқорындағы 17 188 лиганд құрылымымен салыстыру жүргізілуде...",
        "screening_success": "⚡ Скрининг сәтті аяқталды!",
        "top_match_title": "🏆 **Құрылымдық ұқсастық бойынша ең жақсы сәйкестік:** Қапшық / Сайт ID: **{pdb_id}**",
        "top_match_reason": "🧬 *Негіздеме:* Нативті лигандпен Танимото ұқсастық индексі: {sim:.2f}.",
        "top_match_instruction": "💡 *Студентке нұсқаулық:* Сіздің молекулаңыз осы ақуыздың нативті лигандына барынша ұқсас.",
        "table_title": "📊 **Докингке арналған ТОП-15 әлеуетті биологиялық нысаналар:**",
        "col_pdb_id": "PDB ID / Сайт",
        "col_tanimoto": "Танимото индексі (Ұқсастық)",
        "col_pkd": "pKd болжамы",
        "pdb_lookup_title": "🔍 RCSB PDB дерекқорындағы ақуыз құрылымына жылдам өту",
        "pdb_lookup_label": "4 таңбалы PDB ID енгізіңіз (мысалы, 4E5H):",
        "btn_go_to_pdb": "rcsb.org сайтында {pdb_id} нысанасының картасын ашу",
        "pdb_error_length": "⚠️ PDB ID дәл 4 таңбадан тұруы керек.",
        "btn_run_prep": "⚙️ Толық дайындықты бастау",
        "spinner_meeko": "Meeko және RDKit жұмыс істеуде: зарядтар мен торсиондарды есептеу...",
        "docking_success_info": "Лиганд дайын! Торсиондар мен зарядтар есептелді.",
        "docking_error": "PDBQT дайындау кезінде қате кетті.",
        "docking_note_student": "ℹ️ **Студенттерге ескертпе:** Докинг 'кілт пен құлып' принципін имитациялайды. Кілт сәйкес келуі үшін оның байланыс бұрыштары дұрыс болуы керек.",
        "btn_download_pdbqt": "📥 Дайын PDBQT жүктеу",
        "docking_view_label": "🔍 **Дайындалған құрылымды қарау:**",
        "docking_caption": "Оңтайландырылған 3D-модель (талдауға дайындалған).",
        "docking_next_steps_header": "🎓 Әрі қарай не істеу керек?",
        "docking_next_steps_text": """
        Лигандты дайындағаннан кейін сізге:
        1. **Нысана-белокты** дайындау (суды жою, зарядтарды қосу).
        2. **Grid Box** анықтау (белсенді орталықтың координаттары).
        3. AutoDock Vina немесе ұқсас БҚ-да есептеуді іске қосу қажет.
        """,
        "docking_warn_no_3d": "⚠️ Алдымен бірінші қосымша бетте 3D модельді құрастырыңыз!",
               # scPDB target screening / ligand-based target hypotheses
        "docking_main_title": "🧬 Молекулалық докинг: байланысу сайтын таңдау",

        "target_screening_header": "🎯 scPDB бойынша нысана-гипотезалар",
        "target_screening_subtitle": "scPDB нативті лигандтарымен ұқсастық бойынша ықтимал нысаналарды іздеу",
        "target_screening_input_label": "Target-гипотезаларды іздеу үшін SMILES енгізіңіз",
        "target_screening_button": "scPDB бойынша іздеуді бастау",
        "target_screening_top15": "Топ-15 target-гипотеза",
        "target_screening_best_match": "Ең жақсы сәйкестік",
        "target_screening_method": "Әдіс",
        "target_method_short": "scPDB ligand similarity screening",
        "target_method_note": (
            "Әдіс енгізілген молекуланың scPDB базасындағы нативті лигандтармен "
            "ұқсастығына сүйеніп target-гипотезалар ұсынады. Егер кеңейтілген индекс "
            "қолжетімді болса, ranking байланысу сайтының белгілерін, cavity сипаттамаларын "
            "және interaction fingerprint деректерін қосымша ескереді. Нәтиже есептік гипотеза "
            "болып табылады; ол docking есебі, affinity болжамы немесе белсенділіктің "
            "эксперименттік дәлелі емес."
        ),
        "target_error_invalid_smiles": "SMILES қате. Молекула жазбасын тексеріңіз.",
        "target_error_db_unavailable": "scPDB дерекқоры қолжетімсіз немесе бос.",
        "target_error_unknown": "scPDB-скринингті орындау кезінде қате пайда болды.",

        "target_no_hits_message": (
            "scPDB нативті лигандтарымен жеткілікті жақын сәйкестіктер табылмады. "
            "Бұл мишень жоқ дегенді білдірмейді, бірақ ligand-based іздеу күшті гипотеза берген жоқ."
        ),
        "target_candidate_name": "scPDB нысана-гипотезасы: PDB {pdb_id}",
        "target_reason_ligand_similarity": (
            "Енгізілген молекула {pdb_id} кешеніндегі нативті лигандқа ұқсайды. "
            "Танимото индексі = {similarity:.2f}; деңгейі: {similarity_label}."
        ),
        "target_student_interpretation": (
            "Бұл ligand-based гипотеза: егер молекула белок кешенінде табылған лигандқа "
            "ұқсас болса, сәйкес белок әрі қарай талдау үшін ықтимал нысана ретінде қарастырылуы мүмкін."
        ),
        "target_limitation": (
            "Бұл нәтиже docking есебі, affinity болжамы немесе белсенділіктің эксперименттік дәлелі емес."
        ),
        "target_similarity_high": "жоғары ұқсастық",
        "target_similarity_moderate": "орташа ұқсастық",
        "target_similarity_weak": "әлсіз ұқсастық",
        "target_similarity_low": "төмен ұқсастық",
        "target_col_rank": "Ранг",
        "target_col_pdb": "PDB ID",
        "target_col_similarity": "Tanimoto",
        "target_col_similarity_level": "Ұқсастық деңгейі",
        "target_col_score": "Score, %",
        "target_metric_best_similarity": "Ең жақсы Tanimoto ұқсастығы",
        "target_metric_hits": "Порогтан жоғары сәйкестіктер",
        "target_metric_database_size": "scPDB база көлемі",
        "target_min_similarity": "Ұқсастықтың минималды порогы",
        "target_best_match_header": "Ең жақсы сәйкестік",
        "target_screening_stats_header": "scPDB-скрининг қорытындысы",
	"docking_stage3_title": "Docking box және AutoDock Vina конфигурациясы",
	"docking_stage3_desc": "15 scPDB гипотезасының бірін таңдаңыз. Қолданба белсенді орталық координаттарын көрсетеді және AutoDock Vina үшін конфигурация файлдарын дайындайды.",
	"docking_select_candidate_label": "Докинг үшін PDB / scPDB жазбасын таңдаңыз",
	"docking_col_entry_id": "scPDB жазбасы",
	"docking_col_box_center": "Бокс орталығы, Å",
	"docking_col_box_size": "Бокс өлшемі, Å",
	"docking_col_box_status": "Бокс күйі",
	"docking_metric_selected_pdb": "Таңдалған PDB",
	"docking_metric_selected_entry": "scPDB жазбасы",
	"docking_metric_box_status": "Docking box",
	"docking_box_coordinates_title": "#### Docking box координаттары",
	"docking_coord_param": "Параметр",
	"docking_coord_value": "Мәні",
	"docking_box_source_caption": "Орталық мына файлдан есептелді: {center_source}. Өлшем мына файлдан есептелді: {size_source}.",
	"docking_box_unavailable": "Бұл жазба үшін docking box координаттары қолжетімсіз.",
	"docking_config_title": "#### AutoDock Vina үшін дайын docking_config",
	"docking_config_unavailable": "Бұл жазба үшін конфигурация файлы әзірше қолжетімсіз.",
	"btn_download_config_txt": "docking_config.txt жүктеу",
	"btn_download_config_yml": "docking_config.yml жүктеу",
	"docking_receptor_title": "#### Receptor / ақуыз құрылымы",
	"docking_receptor_available": "Бұл жазба үшін дайын receptor.pdbqt қолжетімді.",
	"docking_receptor_not_available": "Дайын receptor.pdbqt қолданбаның бұлттық нұсқасында сақталмайды.",
	"docking_receptor_instruction": "Ақуыз құрылымын RCSB PDB базасынан жүктеп алып, receptor.pdbqt файлын AutoDockTools, Open Babel немесе Meeko арқылы жергілікті компьютерде дайындаңыз. BioSynth-EDU бұлттық нұсқасында тек ықшам scPDB-index, docking box координаттары және конфигурациялар сақталады.",
	"docking_local_command_title": "#### AutoDock Vina жергілікті іске қосу мысалы",
	"docking_requires_molecule": "Алдымен SMILES енгізіп, молекуланы дайындаңыз.",
        "project_warning": "👈 Жобаны бастау үшін бүйірлік панельден қосылысты таңдаңыз (KZ базасы).",
        "project_mol_header": "Зерттеу нысаны",
        "project_smiles_label": "SMILES коды",
        "project_task_header": "Орындалатын тапсырма",
        "project_task_1": "1. Биологиялық белсенділікті болжау (PASS Online).",
        "project_task_2": "2. Липинский ережесіне сәйкестігін бағалау (ADMET қойындысы).",
        "project_task_3": "3. Қорытынды есеп үшін мәліметтерді дайындаңыз.",
        	"project_ketcher_header": "Молекулалық редактор (Ketcher)",
        "project_mod_warning": "Назар аударыңыз: құрылым өзгертілді!",
        "project_btn_analyze": "🧪 Жаңа құрылымды талдау",
        "tab_project": "🚀 Зерттеу жобасы",
        "project_desc": "**Студенттер мен оқушыларға арналған органикалық химиядан зерттеу жобасы**",
        "project_start_hint": "👈 Жұмысты бастау үшін бүйірлік панельден молекуланы таңдаңыз\n(«База kz» бөлімі немесе SMILES қолмен енгізу өрісі)",
        "selected_mol": "Таңдалған қосылыс",
        "input_struct": "🔬 Енгізілген құрылым",
        "current_smiles": "Ағымдағы SMILES",
        "kz_info": "🇰🇿 Қазақстандық әзірлеме туралы мәліметтер",
        "authors": "Авторлар",
        "description": "Сипаттама",
        "classification": "Классификация",
        "year": "Жыл",
        "info_missing": "Ақпарат жоқ",
        "task_title": "📝 Зерттеу жобасына арналған тапсырма",
        "task_full": "Толық тапсырманы ашу",
        "task_step_1": "1. **PASS Online** көмегімен биологиялық белсенділік спектрін болжаңыз.",
        "task_step_2": "2. Фармакологиялық қасиеттерді бағалаңыз (**ADMET** қойындысында) және Липинский ережесін тексеріңіз.",
        "task_step_3": "3. Төмендегі редакторда **молекула құрылымын өзгертіңіз** (функционалдық топтарды қосу/алып тастау, орынбасарларды өзгерту) және бұл өзгерістердің қасиеттерге әсерін 1 және 2-қадамдарды қайталау арқылы талдаңыз.",
        "task_step_4": "4. Баяндама үшін қорытындылар мен ұсыныстар дайындаңыз.",
        "editor_title": "🧪 Молекула құрылымының редакторы",
        "editor_task": "Тапсырма: Төмендегі құрылымды өзгертіп, «Өзгерістерді қолдану» түймесін басыңыз.",
        "structure_changed": "✅ Редакторда құрылым өзгертілді!",
        "new_smiles": "Жаңа SMILES",
        "apply_btn": "🔄 Өзгерістерді қолдану және жобаны жаңарту",
        "change_applied": "✅ Өзгерістер қолданылды!",
        "download_btn": "💾 Өзгертілген SMILES-ті жүктеп алу",
        "change_editor_hint": "👆 Жоғарыдағы редакторда молекуланы өзгертіңіз",
        "lectures_header": "Тақырып бойынша бейне-дәрістер",
        "pass_lecture": "🎥 Дәріс: PASS",
        "adme_lecture": "🎥 Дәріс: ADME",
        "pubchem_lecture": "🎥 Дәріс: PubChem",
        "tools_header": "Талдау құралдары",
        "pass_online": "🌐 PASS Online",
        "swissadme": "🧪 SwissADME",
        "pubchem_search": "📊 PubChem Search",
        "quiz_btn": "📝 Тесттен өтіп, қорғауға арналған сұрақтарды алу",
        "quiz_title": "BioSynth-EDU интеллектуалдық тренажері",
        "quiz_part_1": "1-бөлім: Тестілеу",
        "quiz_part_2": "2-бөлім: Баяндамаға дайындық сұрақтары",
        "check_result": "Нәтижені тексеру",
        "quiz_score": "Сіздің нәтижеңіз",
        "question_label": "Сұрақ",
        "defense_questions": "Баяндамаға дайындалу сұрақтары",
        "close": "Жабу",
        "mol_not_selected": "⚠️ Молекула таңдалмады",
        "how_start": "Қалай бастау керек:",
        "start_step_1": "• Бүйірлік панельдегі Қазақстан каталогынан қосылысты таңдаңыз",
        "start_step_2": "• Немесе SMILES мәнін қолмен енгізіңіз",
        "download_filename": "modified_molecule",
    },
    "English": {
        "nav_title": "BioSynth-EDU",
        "lang_label": "🌐 Select Language",
        "select_placeholder": "-- Select from list --",
        # Sidebar
        "sidebar_select_mol": "🧪 Select Molecule",
        "sidebar_kaz_cat": "🇰🇿 Kazakhstan Developments",
        "sidebar_kaz_label": "Domestic drugs and cases:",
        "sidebar_world_cat": "🌍 Standard Examples",
        "sidebar_world_label": "Drug substance examples:",
        "sidebar_manual": "✍️ Enter SMILES",
        "sidebar_manual_label": "Or paste SMILES below:",
        # Tabs
        "tab_3d": "🔬 Structure & Properties",
        "tab_admet": "📊 ADMET Analysis",
        "tab_docking": "🧬 Docking",
        "tab_project": "🚀 Research Project",
        "tab_edu": "📚 Learning",
        # Categories (World Examples)
        "cat_analgesic": "Analgesic",
        "cat_stimulant": "Stimulant",
        "cat_antipyretic": "Antipyretic",
        "cat_nsaid": "NSAID",
        "cat_antibiotic": "Antibiotic",
        "cat_alkaloid": "Alkaloid",
        "cat_neuro": "Neurotransmitter",
        # Tab 1: 3D & Info
        "btn_build_3d": "🏗️ Build 3D",
        "btn_optimize": "✨ Optimize (MMFF94)",
        "spinner_optimize": "Energy minimization...",
        "style_label": "Display Style",
        "style_surface": "Surface",
        "style_sphere": "Spheres",
        "style_stick": "Sticks",
        "style_line": "Lines",
        "download_sdf": "💾 Download Structure (SDF)",
        "download_help": "SDF file contains 3D atomic coordinates for professional software (AutoDock, PyMOL, etc.)",
        "info_select_mol": "Select a molecule on the left and click 'Build 3D'",
        "header_ref": "📚 Reference",
        "metric_mw": "M. Weight",
        "metric_logp": "LogP",
        "metric_rot_bonds": "Rotatable bonds",
        "spinner_chembl": "Querying ChEMBL...",
        "status_label": "Status:",
        "status_approved": "Approved",
        "mechanism_label": "🔬 Mechanism of Action",
        "no_chembl": "Biological activity not found in ChEMBL",
        "ext_links": "🔗 External Links:",
        "btn_pubchem": "PubChem Profile",
        "btn_chembl": "ChEMBL Data (IC50/Ki)",
        "btn_similarity": "Find Similar",
        "warn_no_pubchem": "Data not found in PubChem",
        "physchem.title": "🔮 Advanced physicochemical analysis",
        "physchem.info.choose_or_enter_structure": "⚠️ Select or enter a molecular structure.",

        "physchem.tab.conformers": "📊 Conformational analysis",
        "physchem.tab.charges": "⚡ Charge distribution",
        "physchem.tab.descriptors": "📝 Topological descriptors",

        "physchem.conformers.description": "**Conformer energy profile analysis (MMFF94; UFF if parameters are unavailable)**",
        "physchem.conformers.num_conformers": "Number of conformers",
        "physchem.conformers.calculate": "Calculate conformational ensemble",
        "physchem.conformers.spinner": "Searching for conformational minima...",
        "physchem.conformers.show_energy_table": "Show energy table",
        "physchem.conformers.download_sdf": "📥 Download minimum-energy 3D structure (SDF)",
        "physchem.conformers.failed": "Could not generate or optimize conformers for this structure. This often happens with invalid SMILES, salts, metal complexes, or molecules without force-field parameters.",
        "physchem.conformers.caption": "Press the button to calculate the conformational ensemble.",
        "physchem.charges.description": "**Partial charge mapping by the Gasteiger-Marsili method**",
        "physchem.charges.calculate": "Calculate and show charges",
        "physchem.charges.spinner": "Calculating 3D geometry and charges...",
        "physchem.charges.heavy_atoms_title": "**Heavy atom charges**",
        "physchem.charges.failed": "Could not calculate charges or build a 3D structure for this molecule.",
        "physchem.charges.caption": "Press the button to calculate charge distribution.",
        "physchem.descriptors.description": "**Reactivity and polar surface descriptors**",
        "physchem.descriptors.failed": "Could not calculate descriptors for this structure.",
        "physchem.table.conformer": "Conformer",
        "physchem.table.conformer_number": "#{number}",
        "physchem.table.energy_kcal_mol": "Energy (kcal/mol)",
        "physchem.table.delta_energy_kcal_mol": "ΔE from minimum (kcal/mol)",
        "physchem.table.atom_index": "Index",
        "physchem.table.atom": "Atom",
        "physchem.table.charge": "Charge",
        "physchem.table.descriptor": "Descriptor",
        "physchem.table.description": "Description",
        "physchem.table.value": "Value",
        "physchem.file.global_minimum_sdf": "global_minimum_structure.sdf",
        "physchem.descriptor.MolWt.name": "MolWt",
        "physchem.descriptor.MolWt.description": "Molecular weight",
        "physchem.descriptor.ExactMolWt.name": "ExactMolWt",
        "physchem.descriptor.ExactMolWt.description": "Exact molecular weight",
        "physchem.descriptor.LogP.name": "LogP",
        "physchem.descriptor.LogP.description": "Crippen logP",
        "physchem.descriptor.MR.name": "MR",
        "physchem.descriptor.MR.description": "Molar refractivity",
        "physchem.descriptor.TPSA.name": "TPSA",
        "physchem.descriptor.TPSA.description": "Topological polar surface area",
        "physchem.descriptor.HBA.name": "HBA",
        "physchem.descriptor.HBA.description": "Hydrogen bond acceptors",
        "physchem.descriptor.HBD.name": "HBD",
        "physchem.descriptor.HBD.description": "Hydrogen bond donors",
        "physchem.descriptor.RotatableBonds.name": "RotatableBonds",
        "physchem.descriptor.RotatableBonds.description": "Rotatable bonds",
        "physchem.descriptor.RingCount.name": "RingCount",
        "physchem.descriptor.RingCount.description": "Ring count",
        "physchem.descriptor.AromaticRings.name": "AromaticRings",
        "physchem.descriptor.AromaticRings.description": "Aromatic rings",
        "physchem.descriptor.HeavyAtoms.name": "HeavyAtoms",
        "physchem.descriptor.HeavyAtoms.description": "Heavy atoms",
        "physchem.descriptor.FractionCSP3.name": "FractionCSP3",
        "physchem.descriptor.FractionCSP3.description": "Fraction of sp3 carbon atoms",
        "physchem.descriptor.QED.name": "QED",
        "physchem.descriptor.QED.description": "QED drug-likeness",
        "admet_header": "📊 ADMET Analysis & Interpretation",
        "admet_instructions": """
        To conduct ADMET analysis:
        1. Click the **"Open ADMETlab 3.0"** button below.
        2. Copy your SMILES and paste it into the input field on the website.
        3. Once the calculation is complete, download the result in **CSV** format.
        4. Upload the file using the Upload button for automatic interpretation.
        """,
        "btn_open_admetlab": "🌐 Open ADMETlab 3.0",
        "btn_open_swissadme": "🧪 Alternative: SwissADME",
        "uploader_label": "📥 Upload ADMETlab CSV results",
        "expander_raw": "👁️ View raw CSV data",
        "subheader_interp": "📝 Brief interpretation of key parameters",
        "metric_lipophilicity": "LogP (Lipophilicity)",
        "metric_bbb": "BBB (Permeability)",
        "metric_tox": "Toxicity / hERG",
        "status_optimal": "✅ Optimal.",
        "status_extreme": "⚠️ Extreme values.",
        "status_bbb_yes": "🧠 Crosses the BBB.",
        "status_bbb_no": "🛡️ Low CNS risk.",
        "status_tox_high": "💔 High risk.",
        "status_tox_low": "❤️ Low risk.",
        "header_lipinski": "🧐 Drug-like Rules Compliance",
        "lipinski_success": "🌟 Molecule fully complies with Lipinski's Rule!",
        "lipinski_warn": "⚠️ Lipinski's Rule violations: ",
        "lipinski_info": "💡 Note: 1 violation is allowed to maintain 'drug-likeness'.",
        "error_interp": "❌ Interpretation error: ",
                # Advanced CNS / BBB / QSAR module
        "advanced_cns_title": "🧠 Advanced CNS / BBB / P-gp investigation",
        "advanced_cns_caption": "A standalone module for advanced QSAR/ADMET analysis of molecules",
        "advanced_cns_body": (
            "This module is designed for a deeper SMILES-based investigation of a molecule. "
            "It provides not only an ADMET prediction, but also an explanation of which molecular "
            "properties affect blood-brain barrier penetration, P-gp efflux risk and potential "
            "CNS exposure."
        ),
        "advanced_cns_features_title": "The module includes:",
        "advanced_cns_features": [
            "RDKit descriptor calculation and BBB profile analysis;",
            "interpretation of LogP, TPSA, MW, HBD/HBA, pKa, charge and P-gp effects;",
            "BBB × P-gp matrix explaining passive permeability and active efflux;",
            "What-if laboratory for educational descriptor modification;",
            "student report generation for a molecule;",
            "additional ML/SHAP analysis of RandomForest models."
        ],
        "advanced_cns_note": (
            "The module is intended for educational and research-oriented QSAR/ADMET analysis. "
            "The results are in silico predictions and do not replace experimental validation."
        ),
        "advanced_cns_button": "Open advanced CNS investigation module",
        "docking_header": "🛠️ Ligand Preparation for Docking",
        "docking_mol_ready": "✅ 3D structure detected and ready for processing.",
        "docking_checklist_title": "**Preparation Checklist:**",
        "docking_checklist_items": """
        1. Adding implicit hydrogens.
        2. 3D conformation generation.
        3. Energy minimization (MMFF94 force field).
        4. **Defining active torsion angles (PDBQT).**
        """,
        "btn_run_prep": "⚙️ Run Full Preparation",
        "spinner_meeko": "Meeko and RDKit at work: calculating charges and torsions...",
        "docking_success_info": "Ligand ready! Torsions and charges calculated.",
        "docking_error": "Error during PDBQT preparation.",
        "docking_note_student": "ℹ️ **Note for students:** Docking simulates a 'lock and key' mechanism. For the key to fit, it must have the correct bond angles.",
        "btn_download_pdbqt": "📥 Download Ready PDBQT",
        "docking_view_label": "🔍 **View Prepared Structure:**",
        "docking_caption": "Optimized 3D model (prepared for analysis).",
        "docking_next_steps_header": "🎓 What's Next?",
        "docking_next_steps_text": """
        After preparing the ligand, you need to:
        1. Prepare the **target protein** (remove water, add charges).
        2. Define the **Grid Box** (active site coordinates).
        3. Run the calculation in AutoDock Vina or similar software.
        """,
        "docking_stage2_title": "🎯 Target Screening via scPDB",
        "docking_stage2_desc": "Based on your ligand's structure, the system performs an instant vector QSAR screening against the scPDB database.",
        "btn_run_screening": "🚀 Run In Silico Target Screening",
        "spinner_screening": "Performing comparison against 17,188 scPDB ligand structures...",
        "screening_success": "⚡ Screening completed successfully!",
        "top_match_title": "🏆 **Top Structural Match:** Folder / Site ID: **{pdb_id}**",
        "top_match_reason": "🧬 *Justification:* Tanimoto similarity index with native ligand: {sim:.2f}.",
        "top_match_instruction": "💡 *Student Instruction:* Your molecule shares high structural similarity with the native ligand of this protein.",
        "table_title": "📊 **TOP-15 Potential Biological Targets for Docking:**",
        "col_pdb_id": "PDB ID / Site",
        "col_tanimoto": "Tanimoto Index (Similarity)",
        "col_pkd": "Predicted pKd",
        "pdb_lookup_title": "🔍 Quick Access to RCSB Protein Data Bank",
        "pdb_lookup_label": "Enter 4-character PDB ID (e.g., 4E5H):",
        "btn_go_to_pdb": "Open target {pdb_id} page on rcsb.org",
        "pdb_error_length": "⚠️ PDB ID must be exactly 4 characters long.",
        "docking_warn_no_3d": "⚠️ Build the 3D model on the first tab first!",
                       # scPDB target screening / ligand-based target hypotheses
        "docking_main_title": "🧬 Molecular docking: Binding-site selection",
        "target_screening_header": "🎯 Target hypotheses from scPDB",
        "target_screening_subtitle": "Search for possible targets based on similarity to native scPDB ligands",
        "target_screening_input_label": "Enter SMILES for target-hypothesis search",
        "target_screening_button": "Run scPDB screening",
        "target_screening_top15": "Top-15 target hypotheses",
        "target_screening_best_match": "Best match",
        "target_screening_method": "Method",
        "target_method_short": "scPDB ligand similarity screening",
        "target_method_note": (
            "The method proposes target hypotheses based on similarity between the query molecule "
            "and native scPDB ligands. When the extended index is available, the ranking also uses "
            "binding-site descriptors, cavity features and interaction fingerprint information. "
            "The result is a computational hypothesis, not docking, not an affinity prediction and "
            "not experimental evidence of activity."
        ),
        "target_error_invalid_smiles": "Invalid SMILES. Please check the molecular notation.",
        "target_error_db_unavailable": "The scPDB target database is unavailable or empty.",
        "target_error_unknown": "An error occurred during scPDB screening.",
        "target_no_hits_message": (
            "No sufficiently close matches to native scPDB ligands were found. "
            "This does not mean that the molecule has no targets, but ligand-based screening "
            "did not produce a strong hypothesis."
        ),
        "target_candidate_name": "scPDB target hypothesis: PDB {pdb_id}",
        "target_reason_ligand_similarity": (
            "The query molecule resembles the native ligand of complex {pdb_id}. "
            "Tanimoto index = {similarity:.2f}; level: {similarity_label}."
        ),
        "target_student_interpretation": (
            "This is a ligand-based hypothesis: if a molecule resembles a ligand found in a protein complex, "
            "the corresponding protein can be considered as a possible target for further analysis."
        ),
        "target_limitation": (
            "This result is not docking, not an affinity prediction and not experimental evidence of activity."
        ),
        "target_similarity_high": "high similarity",
        "target_similarity_moderate": "moderate similarity",
        "target_similarity_weak": "weak similarity",
        "target_similarity_low": "low similarity",
        "target_col_rank": "Rank",
        "target_col_pdb": "PDB ID",
        "target_col_similarity": "Tanimoto",
        "target_col_similarity_level": "Similarity level",
        "target_col_score": "Score, %",
        "target_metric_best_similarity": "Best Tanimoto similarity",
        "target_metric_hits": "Hits above threshold",
        "target_metric_database_size": "scPDB database size",
        "target_min_similarity": "Minimum similarity threshold",
        "target_best_match_header": "Best match",
        "target_screening_stats_header": "scPDB screening summary",
        "project_warning": "👈 Please select a compound in the sidebar (KZ Base) to start the project.",
        "project_mol_header": "Research Object",
        "project_smiles_label": "SMILES code",
        "project_task_header": "Project Task",
        "project_task_1": "1. Perform bioactivity prediction (PASS Online).",
        "project_task_2": "2. Evaluate Lipinski's Rule of Five (ADMET tab).",
        "project_task_3": "3. Prepare data for the final report.",
        "project_ketcher_header": "Molecule Editor (Ketcher)",
        "project_mod_warning": "Warning: Structure has been modified!",
        "project_btn_analyze": "🧪 Analyze New Structure",
	"docking_stage3_title": "Docking box and AutoDock Vina configuration",
	"docking_stage3_desc": "Choose one of the 15 scPDB target hypotheses. The app will show binding-site coordinates and prepare configuration files for AutoDock Vina.",
	"docking_select_candidate_label": "Select a PDB / scPDB entry for docking",
	"docking_col_entry_id": "scPDB entry",
	"docking_col_box_center": "Box center, Å",
	"docking_col_box_size": "Box size, Å",
	"docking_col_box_status": "Box status",
	"docking_metric_selected_pdb": "Selected PDB",
	"docking_metric_selected_entry": "scPDB entry",
	"docking_metric_box_status": "Docking box",
	"docking_box_coordinates_title": "#### Docking box coordinates",
	"docking_coord_param": "Parameter",
	"docking_coord_value": "Value",
	"docking_box_source_caption": "Center calculated from: {center_source}. Size calculated from: {size_source}.",
	"docking_box_unavailable": "Docking box coordinates are unavailable for this entry.",
	"docking_config_title": "#### Ready-to-use AutoDock Vina docking_config",
	"docking_config_unavailable": "The configuration file is not available for this entry yet.",
	"btn_download_config_txt": "Download docking_config.txt",
	"btn_download_config_yml": "Download docking_config.yml",
	"docking_receptor_title": "#### Receptor / protein structure",
	"docking_receptor_available": "A prepared receptor.pdbqt is available for this entry.",
	"docking_receptor_not_available": "A prepared receptor.pdbqt is not stored in the cloud version of the app.",
	"docking_receptor_instruction": "Download the protein structure from RCSB PDB and prepare receptor.pdbqt locally using AutoDockTools, Open Babel, or Meeko. The BioSynth-EDU cloud version stores only the compact scPDB index, docking box coordinates, and configuration files.",
	"docking_local_command_title": "#### Example local AutoDock Vina command",
"docking_requires_molecule": "Enter a SMILES string and prepare the molecule first.",
        "tab_project": "🚀 Research Project",
        "project_desc": "**Organic Chemistry Research Project for Students**",
        "project_start_hint": "👈 Select a molecule in the sidebar to start\n(via 'Catalog kz' or manual SMILES input)",
        "selected_mol": "Selected Compound",
        "input_struct": "🔬 Entered Structure",
        "current_smiles": "Current SMILES",
        "kz_info": "🇰🇿 Details of Kazakhstani Development",
        "authors": "Authors",
        "description": "Description",
        "classification": "Classification",
        "year": "Year",
        "info_missing": "Information not available",
        "task_title": "📝 Research Project Task",
        "task_full": "Open full task",
        "task_step_1": "1. Predict the biological activity spectrum using **PASS Online**.",
        "task_step_2": "2. Evaluate pharmacological properties (in the **ADMET** tab) and check Lipinski's Rule.",
        "task_step_3": "3. **Modify the molecule structure** in the editor below (add/remove functional groups, change substituents) and analyze how these changes affect the properties by repeating steps 1 and 2.",
        "task_step_4": "4. Formulate conclusions and recommendations for the report.",
        "editor_title": "🧪 Molecule Structure Editor",
        "editor_task": "Task: Modify the structure below and click 'Apply Changes'.",
        "structure_changed": "✅ Structure changed in the editor!",
        "new_smiles": "New SMILES",
        "apply_btn": "🔄 Apply Changes and Update Project",
        "change_applied": "✅ Changes applied!",
        "download_btn": "💾 Download Modified SMILES",
        "change_editor_hint": "👆 Modify the molecule in the editor above",
        "lectures_header": "Video Lectures",
        "pass_lecture": "🎥 Lecture: PASS",
        "adme_lecture": "🎥 Lecture: ADME",
        "pubchem_lecture": "🎥 Lecture: PubChem",
        "tools_header": "Analysis Tools",
        "pass_online": "🌐 PASS Online",
        "swissadme": "🧪 SwissADME",
        "pubchem_search": "📊 PubChem Search",
        "quiz_btn": "📝 Take the Quiz & Get Defense Questions",
        "quiz_title": "BioSynth-EDU Intelligent Trainer",
        "quiz_part_1": "Part 1: Quiz",
        "quiz_part_2": "Part 2: Questions for Report Preparation",
        "check_result": "Check result",
        "quiz_score": "Your score",
        "question_label": "Question",
        "defense_questions": "Questions for report preparation",
        "close": "Close",
        "mol_not_selected": "⚠️ No molecule selected",
        "how_start": "How to start:",
        "start_step_1": "• Select a compound from the Kazakhstan catalog in the sidebar",
        "start_step_2": "• Or enter a SMILES string manually in the sidebar",
        "download_filename": "modified_molecule",
    }
}
# ---------------------------------------------------------------------
# scPDB target screening translations
# Добавляется отдельным update-блоком, чтобы не сверять вручную все ключи
# внутри основного словаря LANGUAGES.
# ---------------------------------------------------------------------

TARGET_SCREENING_TRANSLATIONS = {
    "Русский": {
        "docking_main_title": "🧬 Молекулярный докинг: Подбор сайта связывания",

        "target_error_invalid_smiles": "Некорректный SMILES. Проверьте запись молекулы.",
        "target_error_db_unavailable": "База данных scPDB недоступна или пуста.",
        "target_error_unknown": "Ошибка при выполнении scPDB-скрининга.",

        "target_method_note": (
            "Метод сравнивает введённую молекулу с нативными лигандами из scPDB. "
            "В текущей версии используются fingerprint структур ligand.mol2. "
            "Активный сайт напрямую не учитывается, поэтому результат является target-гипотезой, "
            "а не доказательством связывания."
        ),

        "target_no_hits_message": (
            "Не найдено достаточно близких совпадений с нативными лигандами scPDB. "
            "Это не означает отсутствие мишеней, но ligand-based поиск не дал сильной гипотезы."
        ),

        "target_candidate_name": "Гипотеза мишени из scPDB: PDB {pdb_id}",

        "target_reason_ligand_similarity": (
            "Введённая молекула похожа на нативный лиганд комплекса {pdb_id}. "
            "Индекс Танимото = {similarity:.2f}; уровень: {similarity_label}."
        ),

        "target_student_interpretation": (
            "Это ligand-based гипотеза: если молекула похожа на лиганд, найденный "
            "в комплексе с белком, соответствующий белок можно рассматривать как "
            "возможную мишень для дальнейшего анализа."
        ),

        "target_limitation": (
            "Этот результат не является docking-расчётом, прогнозом affinity или "
            "экспериментальным доказательством активности."
        ),

        "target_similarity_high": "высокое сходство",
        "target_similarity_moderate": "умеренное сходство",
        "target_similarity_weak": "слабое сходство",
        "target_similarity_low": "низкое сходство",

        "target_col_rank": "Ранг",
        "target_col_similarity_level": "Уровень сходства",
        "target_col_score": "Score, %",

        "target_best_match_header": "Лучшее совпадение",
        "target_screening_stats_header": "Сводка scPDB-скрининга",
        "target_metric_best_similarity": "Лучшее Tanimoto-сходство",
        "target_metric_hits": "Совпадений выше порога",
        "target_metric_database_size": "Размер базы scPDB",
    },

    "Қазақша": {
        "docking_main_title": "🧬 Молекулалық докинг: байланысу сайтын таңдау",

        "target_error_invalid_smiles": "SMILES қате. Молекула жазбасын тексеріңіз.",
        "target_error_db_unavailable": "scPDB дерекқоры қолжетімсіз немесе бос.",
        "target_error_unknown": "scPDB-скринингті орындау кезінде қате пайда болды.",

        "target_method_note": (
            "Әдіс енгізілген молекуланы scPDB базасындағы нативті лигандтармен салыстырады. "
            "Қазіргі нұсқада ligand.mol2 құрылымдарының fingerprint-тері қолданылады. "
            "Белсенді сайт тікелей есепке алынбайды, сондықтан нәтиже нысана-гипотеза болып табылады, "
            "байланысудың дәлелі емес."
        ),

        "target_no_hits_message": (
            "scPDB нативті лигандтарымен жеткілікті жақын сәйкестіктер табылмады. "
            "Бұл мишень жоқ дегенді білдірмейді, бірақ ligand-based іздеу күшті гипотеза берген жоқ."
        ),

        "target_candidate_name": "scPDB нысана-гипотезасы: PDB {pdb_id}",

        "target_reason_ligand_similarity": (
            "Енгізілген молекула {pdb_id} кешеніндегі нативті лигандқа ұқсайды. "
            "Танимото индексі = {similarity:.2f}; деңгейі: {similarity_label}."
        ),

        "target_student_interpretation": (
            "Бұл ligand-based гипотеза: егер молекула белок кешенінде табылған лигандқа "
            "ұқсас болса, сәйкес белок әрі қарай талдау үшін ықтимал нысана ретінде қарастырылуы мүмкін."
        ),

        "target_limitation": (
            "Бұл нәтиже docking есебі, affinity болжамы немесе белсенділіктің эксперименттік дәлелі емес."
        ),

        "target_similarity_high": "жоғары ұқсастық",
        "target_similarity_moderate": "орташа ұқсастық",
        "target_similarity_weak": "әлсіз ұқсастық",
        "target_similarity_low": "төмен ұқсастық",

        "target_col_rank": "Ранг",
        "target_col_similarity_level": "Ұқсастық деңгейі",
        "target_col_score": "Score, %",

        "target_best_match_header": "Ең жақсы сәйкестік",
        "target_screening_stats_header": "scPDB-скрининг қорытындысы",
        "target_metric_best_similarity": "Ең жақсы Tanimoto ұқсастығы",
        "target_metric_hits": "Порогтан жоғары сәйкестіктер",
        "target_metric_database_size": "scPDB база көлемі",
    },

    "English": {
        "docking_main_title": "🧬 Molecular docking: Binding-site selection",

        "target_error_invalid_smiles": "Invalid SMILES. Please check the molecular notation.",
        "target_error_db_unavailable": "The scPDB target database is unavailable or empty.",
        "target_error_unknown": "An error occurred during scPDB screening.",

        "target_method_note": (
            "The method compares the query molecule with native ligands from scPDB. "
            "In the current version, fingerprints generated from ligand.mol2 structures are used. "
            "The binding site is not directly included, so the result is a target hypothesis, "
            "not proof of binding."
        ),

        "target_no_hits_message": (
            "No sufficiently close matches to native scPDB ligands were found. "
            "This does not mean that the molecule has no targets, but ligand-based screening "
            "did not produce a strong hypothesis."
        ),

        "target_candidate_name": "scPDB target hypothesis: PDB {pdb_id}",

        "target_reason_ligand_similarity": (
            "The query molecule resembles the native ligand of complex {pdb_id}. "
            "Tanimoto index = {similarity:.2f}; level: {similarity_label}."
        ),

        "target_student_interpretation": (
            "This is a ligand-based hypothesis: if a molecule resembles a ligand found in a protein complex, "
            "the corresponding protein can be considered as a possible target for further analysis."
        ),

        "target_limitation": (
            "This result is not docking, not an affinity prediction and not experimental evidence of activity."
        ),

        "target_similarity_high": "high similarity",
        "target_similarity_moderate": "moderate similarity",
        "target_similarity_weak": "weak similarity",
        "target_similarity_low": "low similarity",

        "target_col_rank": "Rank",
        "target_col_similarity_level": "Similarity level",
        "target_col_score": "Score, %",

        "target_best_match_header": "Best match",
        "target_screening_stats_header": "scPDB screening summary",
        "target_metric_best_similarity": "Best Tanimoto similarity",
        "target_metric_hits": "Hits above threshold",
        "target_metric_database_size": "scPDB database size",
    },
}


for _lang_name, _items in TARGET_SCREENING_TRANSLATIONS.items():
    if _lang_name in LANGUAGES:
        LANGUAGES[_lang_name].update(_items)
