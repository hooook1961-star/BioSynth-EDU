import streamlit as st
import json
from pathlib import Path
import openai
import os

def load_tutor_knowledge():
    """Загружает bot_knowledge_new.json и catalog.json с использованием абсолютных путей"""
    data = {"kb": {}, "catalog": []}
    try:
        # Находим абсолютный путь к папке, где лежит сам bot.py (app/core)
        current_dir = Path(__file__).parent.absolute()
        # Папка app — на уровень выше
        app_dir = current_dir.parent
        
        # 1. Инструкции и знания о приложении (ищем в app/data/)
        kb_path = app_dir / "data" / "bot_knowledge_new.json"
        if kb_path.exists():
            with open(kb_path, "r", encoding="utf-8") as f:
                data["kb"] = json.load(f)
        else:
            # Лог для отладки в панели управления Streamlit
            print(f"⚠️ Файл не найден: {kb_path}")
        
        # 2. Каталог молекул (ищем в app/data/)
        cat_path = app_dir / "data" / "catalog.json"
        if cat_path.exists():
            with open(cat_path, "r", encoding="utf-8") as f:
                data["catalog"] = json.load(f)
        else:
            print(f"⚠️ Файл не найден: {cat_path}")
                
        return data
             
    except Exception as e:
        st.error(f"❌ Ошибка загрузки данных тьютора: {e}")
        return data


@st.dialog("🤖 ИИ-Тьютор BioSynth-EDU")
def tutor_dialog():
    """Диалоговое окно с ИИ-тьютором"""
    st.write("**Задайте любой вопрос по казахстанским молекулам** из каталога (Алмакаин, Просидол, Арглабин, Казкаин, Рихлокаин, Пиностробин, Салсолин и др.) или по работе приложения.")
    
    if "tutor_messages" not in st.session_state:
        st.session_state.tutor_messages = []
    
    # Показ истории
    for msg in st.session_state.tutor_messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])
    
    # Ввод вопроса
    if prompt := st.chat_input("Ваш вопрос..."):
        st.session_state.tutor_messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)
        
        with st.spinner("Думаю..."):
            data = load_tutor_knowledge()
            response = ask_ai_tutor(prompt, data)
        
        st.session_state.tutor_messages.append({"role": "assistant", "content": response})
        with st.chat_message("assistant"):
            st.markdown(response)


def ask_ai_tutor(user_query, data):
    try:
        kb = data.get("kb", {})
        catalog = data.get("catalog", [])
        
        # Подготовка списка молекул
        molecule_list = []
        for mol in catalog:
            if isinstance(mol, dict):
                name = mol.get("name", "Unknown")
                classification = mol.get("classification", "")
                molecule_list.append(f"{name} ({classification})")
        
        all_molecules_str = "\n".join(molecule_list) if molecule_list else "Список молекул пуст"

        # Поиск данных по выбранной молекуле
        selected_mol = st.session_state.get('selected_mol_name', 'Не выбрана')
        mol_data_context = ""
        
        if selected_mol and selected_mol not in ["Не выбрана", "None", "", None]:
            search_name = str(selected_mol).strip().lower()
            for mol in catalog:
                if str(mol.get("name", "")).strip().lower() == search_name:
                    mol_data_context = f"\n\n**ДАННЫЕ ПО ВЫБРАННОЙ СТУДЕНТОМ МОЛЕКУЛЕ ({selected_mol}):**\n{json.dumps(mol, ensure_ascii=False, indent=2)}"
                    break

        # Проверка наличия API ключа перед запросом
        if "OPENROUTER_API_KEY" not in st.secrets:
            return "❌ Ошибка: В st.secrets не найден ключ OPENROUTER_API_KEY"

        client = openai.OpenAI(
            base_url="https://openrouter.ai/api/v1",
            api_key=st.secrets["OPENROUTER_API_KEY"],
        )

        # 🔥 ИСПРАВЛЕНИЕ: "openrouter/free" — это не имя модели. 
        # Используем конкретную бесплатную рабочую модель (например, легкую Llama 3)
        response = client.chat.completions.create(
            model="meta-llama/llama-3-8b-instruct:free",
            messages=[
                {"role": "system", "content": f"""Ты — профессиональный ИИ-Тьютор образовательной платформы BioSynth-EDU.
Твоя задача — помогать студентам и отвечать строго на основе предоставленной базы знаний.

КАТАЛОГ КАЗАХСТАНСКИХ МОЛЕКУЛ:
{all_molecules_str}
{mol_data_context}

ЗНАНИЯ О РАБОТЕ И ИНТЕРФЕЙСЕ ПРИЛОЖЕНИЯ:
{json.dumps(kb, ensure_ascii=False, indent=2)}

ПРАВИЛА ОТВЕТА:
1. Отвечай строго по предоставленным данным.
2. Если в каталоге или знаниях приложения нет запрашиваемой информации — честно ответь, что в локальной базе данных этого нет.
3. Не придумывай химические свойства, патенты и несуществующий функционал интерфейса."""},
                {"role": "user", "content": user_query}
            ],
            temperature=0.3,
        )

        return response.choices[0].message.content

    except Exception as e:
        return f"❌ Ошибка отправки запроса в ИИ: {str(e)}"
