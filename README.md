# RabiesMolScreen

**RabiesMolScreen** — это мини‑конвейер для молекулярного докинга, построенный по принципам «чистой» архитектуры.  Код разделён на уровни: абстрактные порты (описание интерфейсов), адаптеры для конкретных внешних утилит (AutoDock Vina, Smina, OpenBabel), сервисы для оркестрации шагов подготовки, докинга и рескоринга, слой IO для чтения/записи данных и валидации контракта и шаблоны отчётов на основе Jinja2.

## Возможности

* Подготовка белков и лигандов к докингу.
* Запуск AutoDock Vina и Smina с сохранением лучших поз.
* Рескоринг результатов и преобразование в Parquet.
* Генерация простых HTML‑отчётов и возможность расширять их визуализациями.
* Чёткий контракт на результаты с версией схемы.
* Гибкая архитектура: новые докинговые движки можно подключать через плагины.

## Установка

Проект можно установить с помощью `pip` или через conda.  В корне репозитория есть файл `requirements.txt` и готовая Conda‑конфигурация `environment.yml`.

### Быстрый вариант (pip)

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Использование conda/mamba

```bash
# создаём окружение
mamba env create -f environment.yml
# активируем окружение
mamba activate rabiesmol
```

## Быстрый старт

Подготовьте каталоги с белками и лигандами, затем выполните последовательно четыре шага:

```bash
# 1. подготовка данных
python -m rabiesmol.cli.main prepare --proteins data/proteins --ligands data/ligands --out data/prepared
# 2. докинг
python -m rabiesmol.cli.main dock --receptor data/prepared/proteins/receptor.pdbqt --ligands data/prepared/ligands --out outputs/run1 --engine vina
# 3. рескоринг
python -m rabiesmol.cli.main rescore-cmd --csv outputs/run1/vina_results.csv --out-parquet outputs/run1/results.parquet
# 4. отчёт
python -m rabiesmol.cli.main report --parquet outputs/run1/results.parquet --out reports/report.html
```

### Параметры CLI

Каждая команда поддерживает опцию `--help` для просмотра всех параметров.

* `prepare` — принимает каталоги с исходными белками и лигандами (`--proteins`, `--ligands`) и каталог для вывода.  В результате создаёт поддиректории `proteins` и `ligands` внутри `--out` и помещает туда подготовленные файлы.
* `dock` — принимает подготовленный файл рецептора (`--receptor`) и каталог с лигандами (`--ligands`), а также каталог для результатов (`--out`) и имя движка (`--engine`, по умолчанию `vina`).  Для каждого лиганда создаётся лог, PDBQT‑файл позы и строка в CSV.
* `rescore-cmd` — выполняет рескоринг CSV с результатами и сохраняет Parquet.  В текущей реализации рескоринг равен исходному score, однако архитектура позволяет подключать внешние алгоритмы.
* `report` — читает Parquet и создаёт простой HTML‑отчёт.  Шаблон можно расширять по вашему вкусу.

## Структура проекта

```
rabiesmol/          # исходный код пакета
├── adapters/       # адаптеры для внешних программ (vina, smina)
├── cli/            # CLI на базе Typer
├── domain/         # описания сущностей и контракта
├── io/             # функции для чтения/записи и валидации
├── ports/          # абстрактные интерфейсы
├── report/         # Jinja2‑шаблоны и генератор HTML
├── services/       # координация операций подготовки, докинга, рескоринга
└── utils/          # утилиты для логирования и вызова внешних команд
```

## Дополнительная информация

* Контракт результатов описан в `rabiesmol/domain/models.py` и `config/results.schema.json`.
* Для подключения новых движков реализуйте подкласс `DockingEngine` и зарегистрируйте его в группе entry‑points `rabiesmol.plugins`.
* Пример тестов и настройка CI находится в папке `tests/` и `.github/workflows/ci.yml`.
