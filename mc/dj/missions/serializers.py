from rest_framework import serializers

from .models import Workflow

class WorkflowSerializer(serializers.ModelSerializer):
    class Meta:
        model = Workflow
        fields = ('uuid', 'serialization', 'status', 'created', 'modified',
                  'mission', 'claimed')
        read_only_fields = ('uuid', 'created', 'modified')
