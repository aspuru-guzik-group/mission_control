from rest_framework import serializers

from .models import Mol

class MolSerializer(serializers.ModelSerializer):
    class Meta:
        model = Mol
        fields = ('uuid', 'created', 'modified', 'cml', 'props')
        read_only_fields = ('uuid', 'created', 'modified',)
