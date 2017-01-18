from django.test import TestCase

from ..models import Mol


class MolTestCase(TestCase):
    def test_has_expected_fields(self):
        kwargs = {
            'cml': 'some cml',
            'props': {'some': 'prop'},
        }
        mol = Mol.objects.create(**kwargs)
        for kwarg, value in kwargs.items():
            self.assertEqual(getattr(mol, kwarg), value)
        expected_attrs = ['uuid', 'created', 'modified']
        for attr in expected_attrs: self.assertTrue(hasattr(mol, attr))
