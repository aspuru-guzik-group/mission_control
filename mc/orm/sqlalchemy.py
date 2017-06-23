import collections
import time
import uuid

from sqlalchemy import Column, MetaData, Table
import sqlalchemy.types as _sa_types

from . import custom_sa_types as _custom_sa_types

def generate_schema():
    metadata = MetaData()
    tables = collections.OrderedDict()
    tables['Mission'] = Table(
        generate_table_name('mission'), metadata,
        generate_key_column(),
        generate_label_column(),
        *generate_timestamp_columns()
    )
    tables['Flow'] = Table(
        generate_table_name('flow'), metadata,
        generate_key_column(),
        generate_label_column(),
        generate_serialization_column(),
        generate_status_column(),
        generate_int_column('num_running_tasks'),
        generate_claimed_column(),
        *generate_timestamp_columns()
    )
    tables['Job'] = Table(
        generate_table_name('job'), metadata,
        generate_key_column(),
        generate_label_column(),
        generate_json_column('job_spec'),
        generate_json_column('data'),
        generate_status_column(),
        generate_claimed_column(),
        *generate_timestamp_columns()
    )
    tables['Queue'] = Table(
        generate_table_name('queue'), metadata,
        generate_key_column(),
        generate_label_column(),
        generate_json_column(column_name='queue_spec'),
        *generate_timestamp_columns()
    )
    tables['Lock'] = Table(
        generate_table_name('lock'), metadata,
        generate_key_column('key'),
        generate_key_column('lockee_key', default=None),
        generate_key_column('locker_key', default=None),
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
def generate_key_column(column_name='key', **kwargs):
    return generate_str_column(
        column_name=column_name, length=36,
        **{'primary_key': True, 'default': str_uuid, **kwargs}
    )

def generate_str_column(column_name=None, length=None, **kwargs):
    assert column_name is not None
    return Column(column_name, _sa_types.String(length=length), **kwargs)

def generate_label_column(column_name='label', **kwargs):
    return Column(column_name, _sa_types.String(length=1024),
                  **{'nullable': True}, **kwargs)

def generate_timestamp_columns():
    return [generate_created_column(), generate_modified_column()]

def generate_created_column():
    return Column('created', _sa_types.Integer(), default=_int_time)

def _int_time(): return int(time.time())

def generate_modified_column():
    return Column('modified', _sa_types.Integer(),
                        default=created_or_int_time, onupdate=_int_time)

def created_or_int_time(context):
    return context.current_parameters.get('created') or _int_time()

def generate_serialization_column(column_name='serialization', **kwargs):
    return Column(column_name, _sa_types.Text(), **{'nullable': True, **kwargs})

def generate_status_column(column_name='status', **kwargs):
    return Column(column_name, _sa_types.String(length=32),
                  **{'default': 'PENDING', **kwargs})

def generate_int_column(column_name=None, **kwargs):
    assert column_name is not None
    return Column(column_name, _sa_types.Integer,
                  **{'nullable': True, **kwargs})

def generate_claimed_column(column_name='claimed', **kwargs):
    return Column(column_name, _sa_types.Boolean(),
                  **{'default': False, **kwargs})

def generate_json_column(column_name='json', **kwargs):
    return Column(column_name, _custom_sa_types.JSONEncodedDict(),
                  **{'nullable': True, **kwargs})
