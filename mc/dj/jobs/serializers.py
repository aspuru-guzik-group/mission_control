from rest_framework import serializers

from .models import Job

class JobSerializer(serializers.ModelSerializer):
    class Meta:
        model = Job
        fields = ('uuid', 'name', 'status', 'created', 'modified',)
        read_only_fields = ('uuid', 'created', 'modified',)
