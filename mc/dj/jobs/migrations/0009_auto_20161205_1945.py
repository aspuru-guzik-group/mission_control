# -*- coding: utf-8 -*-
# Generated by Django 1.10.4 on 2016-12-05 19:45
from __future__ import unicode_literals

import django.contrib.postgres.fields.jsonb
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('jobs', '0008_auto_20161130_1447'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='job',
            name='type',
        ),
        migrations.AddField(
            model_name='job',
            name='spec',
            field=django.contrib.postgres.fields.jsonb.JSONField(default=dict),
        ),
    ]