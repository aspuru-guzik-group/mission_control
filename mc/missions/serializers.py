from rest_framework import serializers

from .models import Flow, Job, Queue

default_json_field_kwargs = {'required': False, 'allow_null': True,
                             'default': dict}

missions_serializers = {}

class FlowSerializer(serializers.ModelSerializer):
    class Meta:
        model = Flow
        fields = ('uuid', 'serialization', 'status',
                  'created', 'modified', 'mission', 'claimed', 'label')
        read_only_fields = ('uuid', 'created', 'modified')
missions_serializers[Flow] = FlowSerializer

class JobSerializer(serializers.ModelSerializer):
    job_spec = serializers.JSONField(**default_json_field_kwargs)
    class Meta:
        model = Job
        fields = ('uuid', 'label', 'status', 'created', 'modified', 'job_spec',
                  'data', 'claimed')
        read_only_fields = ('uuid', 'created', 'modified',)
missions_serializers[Job] = JobSerializer

class QueueSerializer(serializers.ModelSerializer):
    queue_spec = serializers.JSONField(**default_json_field_kwargs)
    class Meta:
        model = Queue
        fields = ['uuid', 'label', 'created', 'modified', 'queue_spec']
        read_only_fields = ['uuid', 'created', 'modified']
missions_serializers[Queue] = QueueSerializer
