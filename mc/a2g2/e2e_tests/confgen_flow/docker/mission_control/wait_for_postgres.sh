#!/bin/bash

set -e

host=$MC_DB_HOST
user=$MC_DB_USER

until PGPASSWORD="$MC_DB_PASSWORD" psql -h "$MC_DB_HOST" -U "$MC_DB_USER" -c '\l'; do
  >&2 echo "Postgres is unavailable - sleeping"
  sleep 1
done

>&2 echo "Postgres is up - executing command"
cmd="$@"
exec $cmd
