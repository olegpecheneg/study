#!/usr/bin/env bash
set -euo pipefail
VENV_NAME="venv"
PYTHON=${PYTHON:-python3}

echo "Creating virtualenv '$VENV_NAME' using $PYTHON"
$PYTHON -m venv "$VENV_NAME"
echo "Activate with: source $VENV_NAME/bin/activate"
echo "Upgrading pip and installing requirements..."
"$VENV_NAME/bin/python" -m pip install --upgrade pip
"$VENV_NAME/bin/python" -m pip install -r requirements.txt
echo "Environment ready."
