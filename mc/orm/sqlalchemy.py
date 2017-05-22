import collections
import datetime
import uuid

from sqlalchemy import Column, MetaData, Table
import sqlalchemy.types as _sa_types

def generate_schema():
    metadata = MetaData()
    tables = collections.OrderedDict()
    tables['mission'] = Table(
        generate_table_name('mission'), metadata,
        generate_uuid_column(),
        generate_label_column(),
        *generate_timestamp_columns()
    )
    tables['flow'] = Table(
        generate_table_name('flow'), metadata,
        generate_uuid_column(),
        generate_label_column(),
        generate_serialization_column(),
        generate_status_column(),
        generate_claimed_column(),
        *generate_timestamp_columns()
    )
    tables['job'] = Table(
        generate_table_name('job'), metadata,
        generate_uuid_column(),
        generate_label_column(),
        generate_serialization_column(),
        generate_status_column(),
        generate_claimed_column(),
        *generate_timestamp_columns()
    )
    tables['queue'] = Table(
        generate_table_name('queue'), metadata,
        generate_uuid_column(),
        generate_label_column(),
        generate_serialization_column(column_name='queue_spec'),
        *generate_timestamp_columns()
    )
    schema = {
        'metadata': metadata,
        'tables': tables,
    }
    return schema

def generate_table_name(table_name=None):
    return '{table_ns}_{table_name}'.format(table_ns='mc_models',
                                            table_name=table_name)

def str_uuid(): return str(uuid.uuid4())
def generate_uuid_column(column_name='uuid', **kwargs):
    return Column(column_name, _sa_types.String(length=16),
                  **{'primary_key': True, 'default': str_uuid, **kwargs})

def generate_label_column(column_name='label', **kwargs):
    return Column(column_name, _sa_types.String(length=1024),
                  **{'nullable': True}, **kwargs)

def generate_timestamp_columns():
    return [generate_created_column(), generate_modified_column()]

def generate_created_column(column_name='created', **kwargs):
    return Column(column_name, _sa_types.DateTime(),
                  **{'default': datetime.datetime.now, **kwargs})

def generate_modified_column(column_name='modified', **kwargs):
    return Column(column_name, _sa_types.DateTime(),
                  **{'onupdate': datetime.datetime.now, **kwargs})

def generate_serialization_column(column_name='serialization', **kwargs):
    return Column(column_name, _sa_types.Text(), **{'nullable': True, **kwargs})

def generate_status_column(column_name='status', **kwargs):
    return Column(column_name, _sa_types.String(length=32),
                  **{'default': 'PENDING', **kwargs})

def generate_claimed_column(column_name='claimed', **kwargs):
    return Column(column_name, _sa_types.Boolean(),
                  **{'default': False, **kwargs})
