#!/bin/bash

set -e

PGCMD="PGPASSWORD=\"$MC_DB_PASSWORD\" psql -h \"$MC_DB_HOST\" -U \"$MC_DB_USER\" -c '\l'"

MAX_ATTEMPTS=5
INTERVAL=3
i=0
CONNECTED=0
while [ $i -le "$MAX_ATTEMPTS" ]; do
  i=$(($i+1))
  eval $PGCMD
  if [ $? -eq 0 ]; then
    CONNECTED=1
    break
  else
    echo "Postgres is unavailable - sleeping"
    sleep $INTERVAL
  fi
done

if [ $CONNECTED -eq 1 ]; then
  echo "Postgres is up - executing command"
  cmd="$@"
  exec $cmd
else
  echo "Failed to connect to postgres, aborting."
  exit 1
fi
