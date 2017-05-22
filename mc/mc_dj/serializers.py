from rest_framework import serializers as _serializers

from mc.orm import sqlalchemy as _mc_sa
from .sa_utils import sa_column_to_serializer_field

sa_schema = _mc_sa.generate_schema()

common_read_only_fields = ['created', 'modified']

flow_table = sa_schema['tables']['flow']
class FlowSerializer(_serializers.Serializer):
    uuid = sa_column_to_serializer_field(flow_table.columns['uuid'])
    label = sa_column_to_serializer_field(flow_table.columns['label'])
    serialization = sa_column_to_serializer_field(
        flow_table.columns['serialization'])
    status = sa_column_to_serializer_field(flow_table.columns['status'])
    claimed = sa_column_to_serializer_field(flow_table.columns['claimed'])
    created = sa_column_to_serializer_field(flow_table.columns['created'])
    modified = sa_column_to_serializer_field(flow_table.columns['modified'])

    class Meta:
        fields = ['uuid', 'label', 'serialization', 'status', 'claimed',
                  'created', 'modified']
        read_only_fields = common_read_only_fields

job_table = sa_schema['tables']['job']
class JobSerializer(_serializers.Serializer):
    uuid = sa_column_to_serializer_field(job_table.columns['uuid'])
    label = sa_column_to_serializer_field(job_table.columns['label'])
    serialization = sa_column_to_serializer_field(job_table.columns['serialization'])
    status = sa_column_to_serializer_field(job_table.columns['status'])
    claimed = sa_column_to_serializer_field(job_table.columns['claimed'])
    created = sa_column_to_serializer_field(job_table.columns['created'])
    modified = sa_column_to_serializer_field(job_table.columns['modified'])

    class Meta:
        fields = ['uuid', 'label', 'serialization', 'status', 'claimed',
                  'created', 'modified']
        read_only_fields = common_read_only_fields

queue_table = sa_schema['tables']['queue']
class QueueSerializer(_serializers.Serializer):
    uuid = sa_column_to_serializer_field(queue_table.columns['uuid'])
    label = sa_column_to_serializer_field(queue_table.columns['label'])
    queue_spec = sa_column_to_serializer_field(
        queue_table.columns['queue_spec'])
    created = sa_column_to_serializer_field(queue_table.columns['created'])
    modified = sa_column_to_serializer_field(queue_table.columns['modified'])

    class Meta:
        fields = ['uuid', 'label', 'queue_spec', 'created', 'modified']
        read_only_fields = common_read_only_fields

def get_model_serializers():
    model_serializers = {}
    from . import models
    class ModelFlowSerializer(_serializers.ModelSerializer, FlowSerializer):
        class Meta(FlowSerializer.Meta):
            model = models.Flow
            read_only_fields = common_read_only_fields
    model_serializers['Flow'] = ModelFlowSerializer

    class ModelJobSerializer(_serializers.ModelSerializer, JobSerializer):
        class Meta(JobSerializer.Meta):
            model = models.Job
    model_serializers['Job'] = ModelJobSerializer

    class ModelQueueSerializer(_serializers.ModelSerializer, QueueSerializer):
        class Meta(QueueSerializer.Meta):
            model = models.Job
    model_serializers['Queue'] = ModelQueueSerializer

    return model_serializers
