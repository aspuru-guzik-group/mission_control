#!/bin/bash

cd $MC_ROOT/mc
./manage.py migrate
./manage.py runserver 0.0.0.0:80
