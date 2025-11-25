Environment setup for DataPreparing

Quick start (Windows PowerShell):

1. Create virtual environment and install requirements:

```powershell
cd d:\pythonProject\MitoFragility\DataPreparing
.\setup_env.ps1
```

2. Activate the environment:

```powershell
.\venv\Scripts\Activate.ps1
```

Quick start (Linux / macOS / cluster with bash):

```bash
cd /path/to/MitoFragility/DataPreparing
./setup_env.sh
source venv/bin/activate
```

Notes:
- Uses `requirements.txt` in this folder.
- Scripts assume the project root is `d:\pythonProject\MitoFragility` or `WORK_DIR` set in `config/pipeline_config.conf`.
