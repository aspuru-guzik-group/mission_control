# -*- coding: utf-8 -*-
# Generated by Django 1.10.3 on 2016-11-30 14:47
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('jobs', '0007_job_type'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='job',
            name='finished',
        ),
        migrations.RemoveField(
            model_name='job',
            name='workflow',
        ),
    ]
