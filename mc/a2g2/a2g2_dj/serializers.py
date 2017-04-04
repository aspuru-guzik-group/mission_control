from rest_framework import serializers

from .models import ChemThing

class ChemThingSerializer(serializers.ModelSerializer):
    class Meta:
        model = ChemThing
        fields = ('uuid', 'created', 'modified', 'types', 'props', 'precursors')
