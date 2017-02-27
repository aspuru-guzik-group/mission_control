import pybel
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_cml_for_conformer(conformer=None):
    sdf_conformer = Chem.MolToMolBlock(conformer.GetOwningMol(),
                                              confId=conformer.GetId())
    return to_cml(input=sdf_conformer, input_format='sdf')

def generate_mols_w_conformers(smiles_list=None):
    smiles_list = smiles_list or ['Cc1ccccc1']
    mols = []
    for smiles in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        mols.append(mol)
    return mols


def write_conformers_to_sd_file(_conformers=None, sd_file_path=None):
    with open(sd_file_path, 'w') as conformer_file:
        writer = Chem.SDWriter(conformer_file)
        for conformer in _conformers:
            writer.write(conformer.GetOwningMol(), confId=conformer.GetId())
        writer.close()

def generate_chemthing_for_conformer(conformer=None, props=None):
    chemthing = {
        'cml': generate_cml_for_conformer(conformer=conformer),
        'props': props
    }
    return chemthing

def from_cml(cml=None, output_format=None):
    return pybel.readstring('cml', cml).write(output_format)

def to_cml(input=None, input_format=None):
    return pybel.readstring(input_format, input).write('cml')
