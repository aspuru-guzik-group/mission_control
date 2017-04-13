from rest_framework import serializers

from .models import Flow, Job


class FlowSerializer(serializers.ModelSerializer):
    class Meta:
        model = Flow
        fields = ('uuid', 'serialization', 'spec', 'status',
                  'created', 'modified', 'mission', 'claimed', 'label')
        read_only_fields = ('uuid', 'created', 'modified')

class JobSerializer(serializers.ModelSerializer):
    class Meta:
        model = Job
        fields = ('uuid', 'name', 'status', 'created', 'modified', 'job_spec',
                  'data', 'error')
        read_only_fields = ('uuid', 'created', 'modified',)
