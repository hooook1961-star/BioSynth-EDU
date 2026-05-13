import streamlit as st
import pandas as pd
import random

def get_assessment_data():
    """Загружает тесты и открытые вопросы из Excel в зависимости от языка."""
    lang = st.session_state.get('lang', 'Русский')
    
    # Маппинг колонок в соответствии с вашими файлами
    lang_map = {
        "Русский": {"q_test": "question_ru", "opt_test": "options_ru", "q_open": "Вопрос (RU)"},
        "Қазақша": {"q_test": "question_kz", "opt_test": "options_kz", "q_open": "Сұрақ (KZ)"},
        "English": {"q_test": "question_en", "opt_test": "options_en", "q_open": "Question (EN)"}
    }
    cols = lang_map.get(lang)

    try:
        # Чтение EXCEL (убедитесь, что пути к файлам верны)
        df_tests = pd.read_excel('data/Тесты.xlsx')
        df_open = pd.read_excel('data/Открытые вопросы.xlsx')

        # Выборка случайных данных: 10 тестов и 3 вопроса к защите
        sampled_tests = df_tests.sample(n=min(10, len(df_tests))).to_dict('records')
        sampled_open = df_open.sample(n=min(3, len(df_open)))[cols['q_open']].tolist()

        return sampled_tests, sampled_open, cols
    except Exception as e:
        st.error(f"Ошибка загрузки Excel: {e}")
        return None, None, None

@st.dialog("BioSynth-EDU Assessment", width="large")
def run_quiz_dialog(t_data, o_qs, c_map, t_dict):
    """Отрисовывает интерфейс тестирования."""
    st.write(f"### {t_dict.get('quiz_part_1', 'Часть 1: Тестирование')}")

    with st.form("quiz_form"):
        user_answers = []

        for i, item in enumerate(t_data):
            q_text = item[c_map['q_test']]
            # Разбиваем варианты (правильный всегда первый в Excel)
            raw_options = [opt.strip() for opt in str(item[c_map['opt_test']]).split(';')]
            correct_answer = raw_options[0]

            # Перемешиваем варианты один раз за сессию диалога
            state_key = f"quiz_shuffle_{hash(q_text)}"
            if state_key not in st.session_state:
                st.session_state[state_key] = random.sample(raw_options, len(raw_options))

            st.write(f"**{i+1}. {q_text}**")
            ans = st.radio(
                f"Radio_{i}",
                options=st.session_state[state_key],
                key=f"radio_ui_{state_key}",
                index=None,
                label_visibility="collapsed"
            )
            user_answers.append((ans, correct_answer))
            st.divider()

        submit_quiz = st.form_submit_button(
            t_dict.get("check_result", "Проверить результат"),
            use_container_width=True
        )

    if submit_quiz:
        # Проверка на то, что все вопросы отвечены
        if any(a is None for a, c in user_answers):
            st.warning("Пожалуйста, ответьте на все вопросы.")
        else:
            score = sum(1 for ans, correct in user_answers if ans == correct)
            st.success(f"{t_dict.get('quiz_score', 'Ваш результат')}: {score} / {len(t_data)}")
            
            st.write(f"### {t_dict.get('quiz_part_2', 'Часть 2: Вопросы к защите')}")
            for q in o_qs:
                st.info(f"❓ {q}")

            if st.button(t_dict.get("close", "Закрыть"), use_container_width=True):
                # Чистим временные ключи перемешивания
                for k in list(st.session_state.keys()):
                    if "quiz_shuffle_" in k:
                        del st.session_state[k]
                st.rerun()
