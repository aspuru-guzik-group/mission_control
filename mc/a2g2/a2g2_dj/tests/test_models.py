from django.test import TestCase

from ..models import ChemThing


class ChemThingTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'types': ['type1', 'type2'],
            'props': {'some': 'prop'},
        }
        chemthing = ChemThing.objects.create(**kwargs)
        for kwarg, value in kwargs.items():
            self.assertEqual(getattr(chemthing, kwarg), value)
        expected_attrs = ['uuid', 'created', 'modified']
        for attr in expected_attrs: self.assertTrue(hasattr(chemthing, attr))
