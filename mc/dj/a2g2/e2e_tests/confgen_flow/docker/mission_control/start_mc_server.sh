#!/bin/bash

cd /mc/dj
./manage.py migrate
./manage.py runserver 0.0.0.0:80
