# app/core/bot.py
import streamlit as st
import json
from pathlib import Path

def load_tutor_knowledge():
    """Загружает bot_knowledge_new.json и catalog.json"""
    data = {"kb": {}, "catalog": []}
    try:
        # 1. Инструкции и знания о приложении
        kb_path = Path("data/bot_knowledge_new.json")
        if kb_path.exists():
            with open(kb_path, "r", encoding="utf-8") as f:
                data["kb"] = json.load(f)
        
        # 2. Каталог молекул
        cat_path = Path("data/catalog.json")
        if cat_path.exists():
            with open(cat_path, "r", encoding="utf-8") as f:
                data["catalog"] = json.load(f)
                
        return data
             
    except Exception as e:
        st.error(f"❌ Ошибка загрузки данных: {e}")
        return data


@st.dialog("🤖 ИИ-Тьютор BioSynth-EDU")
def tutor_dialog():
    """Диалоговое окно с ИИ-тьютором"""
    st.write("Задайте любой вопрос по молекулам, пособию или работе приложения.")
    
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
                    mol_data_context = f"\n\n**ДАННЫЕ ПО {selected_mol}:**\n{json.dumps(mol, ensure_ascii=False, indent=2)}"
                    break

        client = openai.OpenAI(
            base_url="https://openrouter.ai/api/v1",
            api_key=st.secrets["OPENROUTER_API_KEY"],
        )

        response = client.chat.completions.create(
            model="openrouter/free",
            messages=[
                {"role": "system", "content": f"""Ты — ИИ-Тьютор BioSynth-EDU.

КАТАЛОГ МОЛЕКУЛ:
{all_molecules_str}

{mol_data_context}

ЗНАНИЯ О ПРИЛОЖЕНИИ:
{json.dumps(kb, ensure_ascii=False, indent=2)}

Отвечай точно по данным. Не выдумывай."""},
                {"role": "user", "content": user_query}
            ],
            temperature=0.3,
        )

        return response.choices[0].message.content

    except Exception as e:
        return f"❌ Ошибка: {str(e)}"
