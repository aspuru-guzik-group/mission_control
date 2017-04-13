import argparse
import json
import os

from rdkit import Chem
from rdkit.Chem import AllChem

from mc.a2g2.job_modules.confgen import constants as confgen_constants

class ConfgenFixtures(object):
    def __init__(self):
        self.mols = self.generate_mols()
        self.conformers = self.generate_conformers()
        self.confgen_params = self.generate_confgen_params()
    
    def generate_mols(self):
        smiles_list = ['Cc1ccccc1']
        return [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    def generate_conformers(self, num_confs=1):
        mol = self.mols[0]
        AllChem.EmbedMultipleConfs(mol, numConfs=num_confs)
        return mol.GetConformers()

    def generate_confgen_params(self):
        return {'param_%s' % i: 'value_%s' % i for i in range(3)}

    def conformer_to_xyz(self, conformer=None):
        return self.atoms_to_xyz(
            atoms=self.conformer_to_atoms(conformer),
            comment=self.generate_conformer_comment(conformer)
        )

    def atoms_to_xyz(self, atoms=None, comment=None):
        xyz_lines = []
        xyz_lines.append("%s" % len(atoms))
        xyz_lines.append("%s" % (comment or ""))
        for atom in atoms:
            xyz_lines.append("{element} {x:.4f} {y:.4f} {z:.4f}".format(**atom))
        xyz = "\n".join(xyz_lines)
        return xyz

    def conformer_to_atoms(self, conformer=None):
        atoms = []
        for i, atom in enumerate(conformer.GetOwningMol().GetAtoms()):
            pos = conformer.GetAtomPosition(i)
            coords = {coord: float("%.4f" % getattr(pos, coord))
                      for coord in ['x', 'y', 'z']}
            atoms.append({'element': atom.GetSymbol(), **coords})
        return atoms

    def generate_conformer_comment(self, conformer=None):
        return json.dumps({'rdkit_conformer_id': conformer.GetId()})


    @classmethod
    def handle_confgen_command(cls, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument('--outdir')
        parsed_args, extra_args = parser.parse_known_args(args=args)
        outdir = parsed_args.outdir
        xyz_dir = os.path.join(outdir, 'conformers')
        os.makedirs(xyz_dir, exist_ok=True)
        confgen_fixtures = cls()
        for i, conformer in enumerate(confgen_fixtures.conformers):
            xyz_path = os.path.join(xyz_dir, 'conformer_%s.xyz' % i)
            open(xyz_path, 'w').write(
                confgen_fixtures.conformer_to_xyz(conformer))
        confgen_params_outfile_path = os.path.join(
            outdir, confgen_constants.CONFGEN_PARAMS_FILENAME)
        with open(confgen_params_outfile_path, 'w') as f:
            f.write(json.dumps(confgen_fixtures.confgen_params))
