#!/bin/bash

python3 manage.py migrate
echo "Database migration complete."
echo "Starting CoNSEnsX..."
python3 manage.py runserver 0.0.0.0:8000
