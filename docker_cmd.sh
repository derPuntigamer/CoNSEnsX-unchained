#!/bin/bash

#waiting for postgres
until pg_isready --host=consensx-db --username=consensxuser --dbname=consensx
do
  echo "Waiting for PostgreSQL..."
  sleep 1
done

echo "Postgres is ready, running the migrations..."

python3 manage.py makemigrations consensx
python3 manage.py migrate
echo "Database migration complete."
echo "Starting CoNSEnsX..."
python3 manage.py runserver 0.0.0.0:8000
