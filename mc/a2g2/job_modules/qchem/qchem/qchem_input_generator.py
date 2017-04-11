import textwrap


class QChemInputGenerator(object):
    def __init__(self, params=None):
        self.params = params

    def generate_qchem_input(self):
        blocks = []
        for block_name in ['molecule', 'rem']:
            block_handler_name = 'generate_%s_block' % block_name
            block_handler = getattr(self, block_handler_name)
            blocks.append(str(block_handler()))
        qchem_input = "\n".join(blocks)
        return qchem_input

    def generate_molecule_block(self):
        molecule = self.params['molecule']
        return textwrap.dedent(
            '''
            $molecule
            {charge} {multiplicity}
            {atoms_block}
            $end
            '''
        ).format(
            charge=molecule['charge'],
            multiplicity=molecule['multiplicity'],
            atoms_block=self.generate_atoms_block(atoms=molecule['atoms'])
        )

    def generate_atoms_block(self, atoms=None):
        return "\n".join([
            str(self.generate_atom_block(atom=atom))
            for atom in atoms
        ])

    def generate_atom_block(self, atom=None):
        return " ".join(str(atom[attr]) for attr in ['element', 'x', 'y', 'z'])

    def generate_rem_block(self):
        return "\n".join([
            str(self.generate_rem_item_block(rem_item=rem_item))
            for rem_item in self.params.get('rem', {}).items()
        ])

    def generate_rem_item_block(self, rem_item=None):
        block_parts = [rem_item[attr] for attr in ['variable', 'value']]
        if 'comment' in rem_item: block_parts.append(rem_item['comment'])
        return "\t".join(block_parts)

def generate_qchem_input(params=None):
    return QChemInputGenerator(params=params).generate_qchem_input()
