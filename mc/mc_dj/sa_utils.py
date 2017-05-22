from django.db import models as _dj_models
import sqlalchemy.types as _sa_types
from rest_framework import serializers as _serializers


def sa_column_to_dj_field(sa_column=None):
    field_kwargs = {
        'primary_key': sa_column.primary_key,
        'null': sa_column.nullable,
    }
    if sa_column.default:
        field_kwargs['default'] = getattr(sa_column.default.arg, '__wrapped__',
                                          sa_column.default.arg)
    if type(sa_column.type) is _sa_types.String:
        field_cls = _dj_models.CharField
        field_kwargs['max_length'] = sa_column.type.length
    elif type(sa_column.type) is _sa_types.Text:
        field_cls = _dj_models.TextField
    elif type(sa_column.type) is _sa_types.Boolean:
        field_cls = _dj_models.NullBooleanField
    elif type(sa_column.type) is _sa_types.DateTime:
        field_cls = _dj_models.DateTimeField
        if sa_column.default:
            field_kwargs['auto_now_add'] = True
            try: del field_kwargs['default']
            except: pass
        if sa_column.onupdate: field_kwargs['auto_now'] = True
    return field_cls(**field_kwargs)

def sa_column_to_serializer_field(sa_column=None):
    field_kwargs = {}
    if type(sa_column.type) is _sa_types.String:
        field_cls = _serializers.CharField
        field_kwargs['max_length'] = sa_column.type.length
    elif type(sa_column.type) is _sa_types.Text:
        field_cls = _serializers.CharField
    elif type(sa_column.type) is _sa_types.Boolean:
        field_cls = _serializers.NullBooleanField
    elif type(sa_column.type) is _sa_types.DateTime:
        field_cls = _serializers.DateTimeField
    return field_cls(**field_kwargs)
