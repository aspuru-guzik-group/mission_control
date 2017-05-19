import collections
from django.db import models as _dj_models
import sqlalchemy.types as _sa_types

from mc.orm import sqlalchemy as _mc_sa


def generate_models_from_sa_schema(sa_schema=None):
    sa_schema = sa_schema or _mc_sa.generate_schema()
    models = collections.OrderedDict()
    for table_key, sa_table in sa_schema['tables'].items():
        fields = {
            sa_column.name: sa_column_to_dj_field(sa_column)
            for sa_column in sa_table.columns
        }
        meta_cls = type('Meta', (), {
            'app_label': 'mc.mc_dj',
            'db_table': sa_table.name
        })
        model = type(
            table_key.title(), (object,),
            {**fields, 'Meta': meta_cls, '__module__': __name__}
        )
        models[model.__name__] = model
    return models

def sa_column_to_dj_field(sa_column=None):
    field_kwargs = {
        'primary_key': sa_column.primary_key,
        'null': sa_column.nullable,
    }
    if sa_column.default:
        field_kwargs['default'] = getattr(sa_column.default.arg, '__wrapped__',
                                          sa_column.default.arg)
    if isinstance(sa_column.type, _sa_types.String):
        field_cls = _dj_models.CharField
        field_kwargs.update({'max_length': sa_column.type.length})
    elif isinstance(sa_column.type, _sa_types.Text):
        field_cls = _dj_models.TextField
    elif isinstance(sa_column.type, _sa_types.Boolean):
        field_cls = _dj_models.NullBooleanField
    return field_cls(**field_kwargs)

models = generate_models_from_sa_schema()
Mission = models['Mission']
Flow = models['Flow']
Job = models['Job']
Queue = models['Queue']
