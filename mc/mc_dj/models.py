from django.db import models as _dj_models

from .sa_utils import sa_column_to_dj_field

from mc.orm import sqlalchemy as _mc_sa

sa_schema = _mc_sa.generate_schema()

mc_models = []

mission_table = sa_schema['tables']['mission']
class Mission(_dj_models.Model):
    uuid = sa_column_to_dj_field(mission_table.columns['uuid'])
    label = sa_column_to_dj_field(mission_table.columns['label'])
    created = sa_column_to_dj_field(mission_table.columns['created'])
    modified = sa_column_to_dj_field(mission_table.columns['modified'])
    class Meta:
        db_table = mission_table.name
mc_models.append(Mission)

flow_table = sa_schema['tables']['flow']
class Flow(_dj_models.Model):
    uuid = sa_column_to_dj_field(flow_table.columns['uuid'])
    label = sa_column_to_dj_field(flow_table.columns['label'])
    serialization = sa_column_to_dj_field(flow_table.columns['serialization'])
    status = sa_column_to_dj_field(flow_table.columns['status'])
    claimed = sa_column_to_dj_field(flow_table.columns['claimed'])
    created = sa_column_to_dj_field(flow_table.columns['created'])
    modified = sa_column_to_dj_field(flow_table.columns['modified'])
    class Meta:
        db_table = flow_table.name
mc_models.append(Flow)

job_table = sa_schema['tables']['job']
class Job(_dj_models.Model):
    uuid = sa_column_to_dj_field(job_table.columns['uuid'])
    label = sa_column_to_dj_field(job_table.columns['label'])
    serialization = sa_column_to_dj_field(job_table.columns['serialization'])
    status = sa_column_to_dj_field(job_table.columns['status'])
    claimed = sa_column_to_dj_field(job_table.columns['claimed'])
    created = sa_column_to_dj_field(job_table.columns['created'])
    modified = sa_column_to_dj_field(job_table.columns['modified'])
    class Meta:
        db_table = job_table.name
mc_models.append(Job)

queue_table = sa_schema['tables']['queue']
class Queue(_dj_models.Model):
    uuid = sa_column_to_dj_field(queue_table.columns['uuid'])
    label = sa_column_to_dj_field(queue_table.columns['label'])
    queue_spec = sa_column_to_dj_field(queue_table.columns['queue_spec'])
    created = sa_column_to_dj_field(queue_table.columns['created'])
    modified = sa_column_to_dj_field(queue_table.columns['modified'])
    class Meta:
        db_table = queue_table.name
mc_models.append(Queue)
