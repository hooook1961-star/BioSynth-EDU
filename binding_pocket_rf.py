#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

import os
import numpy as np
import deepchem as dc
import joblib  # Библиотека для сохранения модели для Streamlit
from sklearn.ensemble import RandomForestClassifier
from binding_pocket_datasets import load_pdbbind_pockets

# Фиксируем случайные числа для стабильности результатов
np.random.seed(123)

split = "random"
subset = "core"  # Наш лоадер теперь берет Core Set из архива

print("Шаг 1: Загрузка отфильтрованных данных из архива...")
pdbbind_tasks, pdbbind_datasets, transformers = load_pdbbind_pockets(
    split=split, subset=subset)
train_dataset, valid_dataset, test_dataset = pdbbind_datasets

# Задаем метрику оценки (ROC AUC)
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)

# Настраиваем стандартный классификатор Scikit-Learn
sklearn_model = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=123)
model = dc.models.SklearnModel(sklearn_model)

# Запуск обучения
print("Шаг 2: Обучение модели Случайного Леса (Random Forest)...")
model.fit(train_dataset)

# Оценка точности
print("Шаг 3: Проверка качества модели...")
train_scores = model.evaluate(train_dataset, [metric], transformers)
valid_scores = model.evaluate(valid_dataset, [metric], transformers)

print("Метрики на тренировочной выборке (ROC AUC):", train_scores)
print("Метрики на валидационной выборке (ROC AUC):", valid_scores)

# 🔥 НАШЕ ВНЕДРЕНИЕ: Сохранение чистого «мозга» модели для Streamlit
print("Шаг 4: Экспорт модели для веб-интерфейса...")
pure_rf_model = model.model_instance  # Вытаскиваем чистый sklearn объект из обертки

model_filename = "visual_pocket_model.pkl"
joblib.dump(pure_rf_model, model_filename)

print(f"🎉 УСПЕХ! Готовая модель сохранена в файл: {model_filename}")
