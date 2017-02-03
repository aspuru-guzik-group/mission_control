import glob
import os

import pybel


class ConfgenLoadJobEngine(object):
    def parse_job_dir(self, job_dir=None, includes=None, excludes=None):
        paths_to_parse = self.get_paths_to_parse(job_dir=job_dir,
                                                 includes=includes,
                                                 excludes=excludes)
        chemthing_dicts = []
        for path in paths_to_parse: chemthing_dicts.extend(self.parse_path(path))
        return chemthing_dicts

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
        chemthing_dicts = []
        pb_mols = self.file_to_pb_mols(path=path)
        for i, pb_mol in enumerate(pb_mols):
            chemthing_dict = self.pb_mol_to_chemthing_dict(pb_mol)
            chemthing_dict['props']['mol_idx'] = i
            chemthing_dicts.append(chemthing_dict)
        return chemthing_dicts

    def file_to_pb_mols(self, path=None):
        ext = os.path.splitext(path)[1].lstrip('.')
        return pybel.readfile(ext, path)

    def pb_mol_to_chemthing_dict(self, pb_mol=None):
        chemthing_dict = {'cml': pb_mol.write('cml'), 'props': {}}
        return chemthing_dict

