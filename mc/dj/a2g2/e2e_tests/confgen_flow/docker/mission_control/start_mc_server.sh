#!/bin/bash

source activate mission_control_env
cd /mc/dj
./manage.py migrate
./manage.py runserver 0.0.0.0:80
