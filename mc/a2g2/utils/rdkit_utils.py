import os

from rdkit import Chem
from rdkit.Chem import AllChem

from . import cml_utils


def conformer_to_cml(conformer=None):
    sdf_conformer = Chem.MolToMolBlock(conformer.GetOwningMol(),
                                       confId=conformer.GetId())
    return cml_utils.to_cml(input=sdf_conformer, input_format='sdf')

def generate_mols_w_conformers(smiles_list=None):
    smiles_list = smiles_list or ['Cc1ccccc1']
    mols = []
    for smiles in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        mols.append(mol)
    return mols

def write_conformers_to_sd_file(conformers=None, sd_file_path=None):
    with open(sd_file_path, 'w') as conformer_file:
        writer = Chem.SDWriter(conformer_file)
        for conformer in conformers:
            writer.write(conformer.GetOwningMol(), confId=conformer.GetId())
        writer.close()

def write_conformers_to_xyz_files(conformers=None, parent_dir=None,
                                  xyz_name_tpl=None, extra_tpl_kwargs=None):
    xyz_name_tpl = xyz_name_tpl or 'conformer__rdk_{rdk_id}'
    extra_tpl_kwargs = extra_tpl_kwargs or {}
    xyz_paths = []
    for idx, conformer in enumerate(conformers):
        xyz_str = conformer_to_xyz(conformer)
        tpl_kwargs = {'rdk_id': conformer.GetId(), **extra_tpl_kwargs}
        xyz_name = xyz_name_tpl.format(**tpl_kwargs)
        xyz_path = os.path.join(parent_dir, xyz_name)
        with open(xyz_path, 'w') as xyz_file: xyz_file.write(xyz_str)
        xyz_paths.append(xyz_path)

def conformer_to_xyz(conformer=None, comment=""):
    atoms = conformer_to_atoms(conformer=conformer)
    return atoms_to_xyz_str(atoms=atoms, comment=comment)

def conformer_to_atoms(conformer=None):
    atoms = []
    for i, atom in enumerate(conformer.GetOwningMol().GetAtoms()):
        pos = conformer.GetAtomPosition(i)
        atoms.append([atom.GetSymbol(), [pos.x, pos.y, pos.z]])
    return atoms

def atoms_to_xyz_str(atoms=None, comment=None):
    xyz_lines = []
    xyz_lines.append("%s" % len(atoms))
    xyz_lines.append("%s" % (comment or ""))
    for atom in atoms:
        xyz_lines.append("%s %.4f %.4f %.4f" % (atom[0], *atom[1]))
    xyz_str = "\n".join(xyz_lines)
    return xyz_str
