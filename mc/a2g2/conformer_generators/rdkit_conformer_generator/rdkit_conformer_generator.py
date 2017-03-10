import os
import subprocess
import re
import shutil
import tempfile

from rdkit.Chem.AllChem import (
    AddHs, EmbedMultipleConfs, GetBestRMS,
    MMFFGetMoleculeForceField, MMFFGetMoleculeProperties,
    MolFromSmiles,RemoveHs,
    UFFGetMoleculeForceField
)

class RDKitConformerGenerator(object):

    forcefields = ['mmff', 'uff']

    '''
    Generates low-resolution conformers.
    '''
    def __init__(
        self,
        output_limit=20,
        candidate_limit=200,
        forcefield_id='mmff',
        output_dir=None,
        smiles=None,
        prune_rms_thresh=0.01,
        cluster_energy_range=5.0,
        cluster_energy_radius=5.0,
        cluster_energy_auto_accept_radius=1e-3,
        cluster_rms_radius=0.1,
        xyz_filename_tpl='Conf_{index}.xyz',
        fallback_to_align_for_rms=False,
        **kwargs
    ):
        if 'B' in smiles and forcefield_id == 'mmff':
            raise Exception("Smiles '{}' contains boron, but MMFF94"
                            " does not have boron\n".format(smiles))
        if forcefield_id != "mmff" and forcefield_id != "uff":
            raise Exception(
                "Unrecognised force field '{}'".format(forcefield_id))

        self.output_limit = output_limit
        self.candidate_limit = candidate_limit
        self.smiles = smiles
        self.output_dir = output_dir
        self.prune_rms_thresh = prune_rms_thresh
        self.cluster_energy_range = cluster_energy_range
        self.cluster_energy_radius = cluster_energy_radius
        self.cluster_energy_auto_accept_radius = (
            cluster_energy_auto_accept_radius)
        self.cluster_rms_radius = cluster_rms_radius
        self.forcefield_id = forcefield_id
        self.fallback_to_align_for_rms = fallback_to_align_for_rms
        self.mol = AddHs(MolFromSmiles(self.smiles), addCoords=True)
        self.xyz_filename_tpl = xyz_filename_tpl

    def generate_conformers(self):
        conformer_ids = self._get_conformer_ids()
        energys = self._get_energys(conformer_ids=conformer_ids)
        squashed_energys = self._squash_energys(energys)
        conformers_for_squashed_energys = [
            self.mol.GetConformer(energy['conformer_id'])
            for energy in squashed_energys
        ]
        tmp_dir = tempfile.mkdtemp(prefix="confgen.")
        xyz_files = []
        for conformer in conformers_for_squashed_energys:
            xyz_basename = "conf_id_%s.xyz" % conformer.GetId()
            xyz_file = os.path.join(tmp_dir, xyz_basename)
            with open(xyz_file, 'w') as f:
                f.write(self._conformer_to_xyz(conformer))
            xyz_files.append(xyz_file)
        if not xyz_files:
            raise Exception("No candidate conformers after energy clustering.")

        deannotated_smarts = self._smiles_to_deannotated_smarts(self.smiles)

        rms_squashed_xyz_files = self._squash_xyz_files_by_rms(
            xyz_files=xyz_files,
            smarts=deannotated_smarts
        )
        for i, squashed_xyz_file in enumerate(rms_squashed_xyz_files):
            xyz_filename = self.xyz_filename_tpl.format(index=(i + 1))
            output_xyz_path = os.path.join(self.output_dir, xyz_filename)
            shutil.copy(squashed_xyz_file, output_xyz_path)

    def _get_conformer_ids(self):
        conformer_ids = EmbedMultipleConfs(
            self.mol,
            numConfs=self.candidate_limit,
            pruneRmsThresh=self.prune_rms_thresh,
            useRandomCoords=False,
            randomSeed=1,
        )
        return [int(_id) for _id in conformer_ids]

    def _get_energys(self, conformer_ids=None):
        energys = []
        if self.forcefield_id == "mmff":
            mol_props = MMFFGetMoleculeProperties(self.mol)
            def _get_forcefield(confId):
                return MMFFGetMoleculeForceField(self.mol, mol_props,
                                                 confId=confId)
        elif self.forcefield_id == "uff":
            def _get_forcefield(confId):
                return UFFGetMoleculeForceField(self.mol, confId=confId)
        for conformer_id in conformer_ids:
            forcefield = _get_forcefield(confId=conformer_id)
            forcefield.Minimize()
            energys.append({'conformer_id': conformer_id,
                            'value': forcefield.CalcEnergy()})
        return energys

    def _squash_energys(self, energys=None):
        if not energys: return []

        def _get_energy_delta(energy_a, energy_b):
            return abs(energy_a['value'] - energy_b['value'])

        clustered_indexs = {}
        def _add_to_cluster(cluster, energy, index, rms=None):
            cluster['members'].append({'energy': energy, 'rm': rms })
            clustered_indexs[index] = index

        def _generate_new_cluster(seed, index):
            cluster = {'seed': seed, 'members': []}
            _add_to_cluster(cluster, seed, index, rms=0)
            return cluster

        clusters = []
        sorted_energys = sorted(energys, key=lambda d: d['value'])
        min_energy = sorted_energys[0]
        mol_no_hs = RemoveHs(self.mol)
        atom_map = [(i, i) for i in range(mol_no_hs.GetNumAtoms())]
        num_energys = len(sorted_energys)
        for i in range(num_energys):
            if i in clustered_indexs: continue
            seed = sorted_energys[i]
            if _get_energy_delta(seed, min_energy) > self.cluster_energy_range:
                break
            cluster = _generate_new_cluster(seed, i)
            for j in range(num_energys):
                if j in clustered_indexs: continue
                candidate = sorted_energys[j]
                candidate_delta = _get_energy_delta(seed, candidate)
                if candidate_delta <= self.cluster_energy_auto_accept_radius:
                    _add_to_cluster(cluster, candidate, j)
                    continue
                if candidate_delta <= self.cluster_energy_radius:
                    rms = GetBestRMS(
                        mol_no_hs,
                        mol_no_hs,
                        refConfId=seed['conformer_id'],
                        probeConfId=candidate['conformer_id'],
                        maps=[atom_map]
                    )
                    if rms <= self.cluster_rms_radius:
                        _add_to_cluster(cluster, candidate, j, rms=rms)
                        continue
            clusters.append(cluster)
        squashed_energys = [cluster['seed'] for cluster in clusters]
        return squashed_energys

    def _squash_xyz_files_by_rms(self, xyz_files=None, smarts=None):
        files_to_exclude = {}
        inverted_files = {}
        inverted_files_dir = tempfile.mkdtemp(prefix='confgen.inverted.')
        num_files = len(xyz_files)
        for i in range(num_files):
            xyz_a_file = xyz_files[i]
            if xyz_a_file in files_to_exclude: continue
            for j in range(i + 1, num_files):
                xyz_b_file = xyz_files[j]
                if xyz_b_file in files_to_exclude: continue
                rms = self._get_rms(
                    xyz_a_file=xyz_a_file,
                    xyz_b_file=xyz_b_file,
                    smarts=smarts,
                )
                if rms > self.cluster_rms_radius:
                    if xyz_b_file not in inverted_files:
                        with open(xyz_b_file) as f:
                            xyz_b = f.read()
                        inverted_xyz_b = self._invert_xyz(xyz_b)
                        inverted_xyz_b_file = os.path.join(
                            inverted_files_dir, os.path.basename(xyz_b_file))
                        with open(inverted_xyz_b_file, 'w') as f:
                            f.write(inverted_xyz_b)
                        inverted_files[xyz_b_file] = inverted_xyz_b_file
                    inverted_xyz_b_file = inverted_files[xyz_b_file]
                    rms = self._get_rms(
                        xyz_a_file=xyz_a_file,
                        xyz_b_file=inverted_xyz_b_file,
                        smarts=smarts,
                    )
                if rms <= self.cluster_rms_radius:
                    files_to_exclude[xyz_b_file] = xyz_b_file
        squashed_xyz_files = [f for f in xyz_files
                              if f not in files_to_exclude]
        return squashed_xyz_files

    def _conformer_to_xyz(self, conformer, comment=""):
        atoms = self._conformer_to_atoms(conformer)
        return self._atoms_to_xyz_str(atoms=atoms, comment=comment)

    def _conformer_to_atoms(self, conformer):
        atoms = []
        for i, atom in enumerate(conformer.GetOwningMol().GetAtoms()):
            pos = conformer.GetAtomPosition(i)
            atoms.append([atom.GetSymbol(), [pos.x, pos.y, pos.z]])
        return atoms

    def _atoms_to_xyz_str(self, atoms=None, comment=None):
        xyz_lines = []
        xyz_lines.append("%s" % len(atoms))
        xyz_lines.append("%s" % (comment or ""))
        for atom in atoms:
            xyz_lines.append("%s %.4f %.4f %.4f" % (atom[0], *atom[1]))
        xyz_str = "\n".join(xyz_lines)
        return xyz_str

    def _invert_xyz(self, xyz):
        xyz_lines = xyz.split("\n")
        comment = xyz_lines[1]
        inverted_atoms = []
        for atom in self._xyz_str_to_atoms(xyz):
            inverted_atom = [atom[0], [-1 * coord for coord in atom[1]]]
            inverted_atoms.append(inverted_atom)
        return self._atoms_to_xyz_str(inverted_atoms, comment=comment)

    def _xyz_str_to_atoms(self, xyz_str):
        atoms = []
        for atom_line in xyz_str.split("\n")[2:]:
            atom_parts = atom_line.split()
            element = atom_parts[0]
            coords = [float(coord) for coord in atom_parts[1:]]
            atoms.append([element, coords])
        return atoms

    def _xyz_file_to_smiles(self, xyz_file):
        cmd = ["obabel", xyz_file, "-osmi"]
        smiles = self._run_command(cmd).split()[0]
        return smiles

    def _smiles_to_deannotated_smarts(self, smiles):
        deannotated_smarts = smiles
        subs = [
            ('C', '[#6]',),
            ('c', '[#6]',),
            ('N', '[#7]',),
            ('n', '[#7]',),
            ('\[\[', '[',),
            ('\]\]', ']',),
            ('\]H\]', 'H]',),
            ('=', '~',),
        ]
        for sub in subs:
            deannotated_smarts = re.sub(*sub, deannotated_smarts)
        return deannotated_smarts

    def _get_rms(self, xyz_a_file=None, xyz_b_file=None, smarts=None):
        try:
            rms = self._get_obfit_rms(xyz_a_file, xyz_b_file, smarts)
        except subprocess.CalledProcessError:
            if self.fallback_to_align_for_rms:
                rms = self._get_align_rms(xyz_a_file, xyz_b_file)
            else:
                raise
        return rms

    def _get_obfit_rms(self, xyz_a_file, xyz_b_file, smarts):
        cmd = ["obfit", smarts, xyz_a_file, xyz_b_file]
        result = self._run_command(cmd)
        rms = float(result[5:13])
        return rms

    def _get_align_rms(self, xyz_a_file, xyz_b_file, smarts=None):
        cmd = ["obabel", xyz_a_file, xyz_b_file,
               '-o', 'smi', '--align', '--append', 'rmsd']
        if smarts:
            cmd += ['-s', smarts]
        result = self._run_command(cmd)
        rms = float(result.split()[-4])
        return rms

    def _run_command(self, command):
        return subprocess.check_output(
            command, stdin=None,
            stderr=subprocess.STDOUT,
            shell=False,
            universal_newlines=False).decode('utf-8')
