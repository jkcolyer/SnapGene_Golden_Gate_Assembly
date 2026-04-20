#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON=python3

# Check if streamlit is available in this Python
if ! $PYTHON -m streamlit --version &>/dev/null; then
    echo "Streamlit not found in this Python environment. Installing..."
    $PYTHON -m pip install -r requirements.txt
fi

echo "Starting Golden Gate Assembly app..."
echo "Open your browser to http://localhost:8501"

python3 -m streamlit run app.py