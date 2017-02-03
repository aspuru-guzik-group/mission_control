import glob
import os

import pybel


class ConfgenLoadJobEngine(object):
    def __init__(self, a2g2_client=None):
        self.a2g2_client = a2g2_client

    def execute_job(self, job=None):
        chemthings = self.parse_job_dir(job_dir=job['input']['dir_to_parse'])
        job['output'] = self.upload_chemthings(chemthings)

    def parse_job_dir(self, job_dir=None, includes=None, excludes=None):
        paths_to_parse = self.get_paths_to_parse(job_dir=job_dir,
                                                 includes=includes,
                                                 excludes=excludes)
        chemthings = []
        for path in paths_to_parse: chemthings.extend(self.parse_path(path))
        return chemthings

    def get_paths_to_parse(self, job_dir=None, includes=None, excludes=None):
        if not includes: includes = ['raw_data/**']
        paths_to_include = self.get_paths_for_globs(root_dir=job_dir,
                                                    globs=includes)
        paths_to_exclude = self.get_paths_for_globs(root_dir=job_dir,
                                                    globs=excludes)
        paths_to_parse = set(paths_to_include).difference(set(paths_to_exclude))
        return paths_to_parse

    def get_paths_for_globs(self, root_dir=None, globs=None):
        paths = []
        globs = globs or []
        for _glob in globs:
            full_glob = os.path.join(root_dir, _glob)
            all_paths = glob.glob(full_glob, recursive=True)
            file_paths = [path for path in all_paths if os.path.isfile(path)]
            paths.extend(file_paths)
        return paths

    def parse_path(self, path=None):
        chemthings = []
        pb_mols = self.file_to_pb_mols(path=path)
        for i, pb_mol in enumerate(pb_mols):
            chemthing = self.pb_mol_to_chemthing(pb_mol)
            chemthing['props']['mol_idx'] = i
            chemthings.append(chemthing)
        return chemthings

    def file_to_pb_mols(self, path=None):
        ext = os.path.splitext(path)[1].lstrip('.')
        return pybel.readfile(ext, path)

    def pb_mol_to_chemthing(self, pb_mol=None):
        chemthing = {'cml': pb_mol.write('cml'), 'props': {}}
        return chemthing

    def upload_chemthings(self, chemthings=None):
        results = []
        for chemthing in chemthings:
            result = self.a2g2_client.create_chemthing(chemthing=chemthing)
            results.append(result)
        return results
