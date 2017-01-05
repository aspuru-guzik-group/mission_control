from rest_framework import serializers

from .models import Flow

class FlowSerializer(serializers.ModelSerializer):
    class Meta:
        model = Flow
        fields = ('uuid', 'serialization', 'status', 'created', 'modified',
                  'mission', 'claimed')
        read_only_fields = ('uuid', 'created', 'modified')
