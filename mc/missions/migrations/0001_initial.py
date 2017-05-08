# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-05-05 12:56
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion
import django_extensions.db.fields
import jsonfield.fields
import missions.models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Flow',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False)),
                ('label', models.CharField(blank=True, max_length=256, null=True)),
                ('serialization', models.TextField(null=True)),
                ('spec', models.TextField(null=True)),
                ('status', models.CharField(choices=[('PENDING', 'pending'), ('RUNNING', 'running'), ('COMPLETED', 'completed'), ('FAILED', 'failed')], default='PENDING', max_length=32, null=True)),
                ('claimed', models.NullBooleanField(default=False)),
            ],
            options={
                'get_latest_by': 'modified',
                'ordering': ('-modified', '-created'),
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='FlowJob',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False)),
                ('finished', models.NullBooleanField(default=False)),
                ('meta', jsonfield.fields.JSONField(default=dict)),
                ('flow', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='flow_jobs', to='missions.Flow')),
            ],
            options={
                'get_latest_by': 'modified',
                'ordering': ('-modified', '-created'),
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Job',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=1024, null=True)),
                ('job_spec', jsonfield.fields.JSONField(default=dict)),
                ('status', models.CharField(choices=[('PENDING', 'pending'), ('RUNNING', 'running'), ('COMPLETED', 'completed'), ('FAILED', 'failed')], default='PENDING', max_length=32, null=True)),
                ('claimed', models.NullBooleanField(default=False)),
                ('data', models.TextField(default='{}', null=True)),
                ('error', models.TextField(blank=True, null=True)),
            ],
            options={
                'get_latest_by': 'modified',
                'ordering': ('-modified', '-created'),
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Mission',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=1024, null=True)),
            ],
            options={
                'get_latest_by': 'modified',
                'ordering': ('-modified', '-created'),
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Queue',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.CharField(default=missions.models.str_uuid4, editable=False, max_length=64, primary_key=True, serialize=False)),
                ('label', models.CharField(blank=True, max_length=256, null=True)),
                ('queue_spec', jsonfield.fields.JSONField(blank=True, default=dict, null=True)),
            ],
            options={
                'get_latest_by': 'modified',
                'ordering': ('-modified', '-created'),
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='flowjob',
            name='job',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='missions.Job'),
        ),
        migrations.AddField(
            model_name='flow',
            name='mission',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='missions.Mission'),
        ),
    ]
