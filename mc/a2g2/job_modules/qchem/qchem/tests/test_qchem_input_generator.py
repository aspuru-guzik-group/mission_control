import collections
import textwrap
import unittest
from unittest.mock import call, MagicMock

from .. import qchem_input_generator

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.params = collections.defaultdict(MagicMock)
        self.generator = self.generate_generator(params=self.params)

    def generate_generator(self, params=None):
        return qchem_input_generator.QChemInputGenerator(params=params)

class GenerateQChemInputTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mock_generator_methods()

    def mock_generator_methods(self):
        for block_name in ['molecule', 'rem']:
            setattr(self.generator, 'generate_%s_block' % block_name,
                    MagicMock())

    def test_generates_expected_input(self):
        blocks = []
        for block_name in ['molecule', 'rem']:
            block_handler_name = 'generate_%s_block' % block_name
            block_handler = getattr(self.generator, block_handler_name)
            blocks.append(str(block_handler()))
        expected_qchem_input = "\n".join(blocks)
        self.assertEqual(self.generator.generate_qchem_input(),
                         expected_qchem_input)

class GenerateMoleculeBockTestCase(BaseTestCase):
    def test_generates_expected_molecule_block(self):
        self.generator.generate_atoms_block = MagicMock()
        expected_molecule_block = textwrap.dedent(
            '''
            $molecule
            {charge} {multiplicity}
            {atoms_block}
            $end
            '''
        ).format(
            charge=self.generator.params['molecule']['charge'],
            multiplicity=self.generator.params['molecule']['multiplicity'],
            atoms_block=self.generator.generate_atoms_block.return_value,
        )
        self.assertEqual(self.generator.generate_molecule_block(),
                         expected_molecule_block)
        self.assertEqual(self.generator.generate_atoms_block.call_args,
                         call(atoms=self.generator.params['molecule']['atoms']))

class GenerateAtomsBlockTestCase(BaseTestCase):
    def test_generates_expected_atoms_block(self):
        self.generator.generate_atom_block = MagicMock()
        atoms = [MagicMock for i in range(3)]
        expected_atoms_block = "\n".join([
            str(self.generator.generate_atom_block(atom=atom))
            for atom in atoms
        ])
        self.assertEqual(self.generator.generate_atoms_block(atoms=atoms),
                         expected_atoms_block)

class GenerateAtomBlockTestCase(BaseTestCase):
    def test_generates_expected_atom_line(self):
        atom = {attr: '%s_value' % attr for attr in ['element', 'x', 'y', 'z']}
        expected_atom_block = " ".join(atom[attr]
                                       for attr in ['element', 'x', 'y', 'z'])
        self.assertEqual(self.generator.generate_atom_block(atom=atom),
                         expected_atom_block)

class GenerateRemBlockTestCase(BaseTestCase):
    def test_generates_expected_rem_block(self):
        self.generator.generate_rem_item_block = MagicMock()
        self.generator.params['rem'] = {'rem_item_%s' % i: MagicMock()
                                        for i in range(3)}
        expected_rem_block = "\n".join([
            str(self.generator.generate_rem_item_block(rem_item=rem_item))
            for rem_item in self.generator.params['rem'].items()
        ])
        self.assertEqual(self.generator.generate_rem_block(),
                         expected_rem_block)

class GenerateRemItemBlockTestCase(BaseTestCase):
    def test_generates_expected_rem_item_block(self):
        rem_item_sans_comment = {'variable': 'some_variable',
                                 'value': 'some_value'}
        rem_item_with_comment = {**rem_item_sans_comment,
                                 'comment': 'some_comment'}
        rem_items = [rem_item_sans_comment, rem_item_with_comment]
        actual_blocks = []
        expected_blocks = []
        for rem_item in rem_items:
            actual_blocks.append(self.generator.generate_rem_item_block(
                rem_item=rem_item))
            expected_block_parts = [rem_item['variable'], rem_item['value']]
            if 'comment' in rem_item:
                expected_block_parts.append(rem_item['comment'])
            expected_blocks.append("\t".join(expected_block_parts))
        self.assertEqual(actual_blocks, expected_blocks)
