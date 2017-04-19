from rest_framework import serializers

from .models import ChemThing

default_json_field_kwargs = {'required': False, 'allow_null': True,
                             'default': dict}

class ChemThingSerializer(serializers.ModelSerializer):
    types = serializers.JSONField(**default_json_field_kwargs)
    keys = serializers.JSONField(**default_json_field_kwargs)
    props = serializers.JSONField(**default_json_field_kwargs)
    precursors = serializers.JSONField(**default_json_field_kwargs)
    ancestors = serializers.JSONField(**default_json_field_kwargs)

    class Meta:
        model = ChemThing
        fields = ('uuid', 'created', 'modified', 'types', 'props', 'precursors',
                  'ancestors', 'keys')
