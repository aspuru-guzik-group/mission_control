# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-01-12 21:48
from __future__ import unicode_literals

from django.db import migrations, models
import missions.models


class Migration(migrations.Migration):

    dependencies = [
        ('missions', '0017_flow_claimed'),
    ]

    operations = [
        migrations.AddField(
            model_name='flow',
            name='spec',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='flow',
            name='status',
            field=models.CharField(choices=[('PENDING', 'pending'), ('RUNNING', 'running'), ('COMPLETED', 'completed')], default='PENDING', max_length=32, null=True),
        ),
        migrations.AlterField(
            model_name='flow',
            name='uuid',
            field=models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='flowjob',
            name='uuid',
            field=models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='mission',
            name='uuid',
            field=models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False),
        ),
    ]