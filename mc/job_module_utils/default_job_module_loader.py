import importlib
import random

import dill


class DefaultJobModuleLoader(object):
    """A class that helps with loading job modules."""

    class JobModuleImportError(Exception): pass

    def __init__(self):
        self.overrides = {}

    def load_job_module(self, job=None, cfg=None):
        """
        Args:
            job (dict): job dict.
            cfg (dict): cfg dict.

        Returns:
            job_module (module): the loaded job module.
        """
        try:
            job_module_name = self.get_job_module_name(job=job, cfg=cfg)
            override = self.overrides.get(job_module_name)
            if override:
                module = self.load_module_per_override(
                    job=job, cfg=cfg, override=override)
            else:
                module = self.load_module_from_module_path(
                    module_path=job_module_name)
            return module
        except Exception as exc:
            msg = "Could not load module for job"
            raise self.JobModuleImportError(msg) from exc

    def load_module_per_override(self, job=None, cfg=None, override=None):
        if isinstance(override, str):
            override = {'type': 'py_module_path', 'params': {'path': override}}
        if override['type'] == 'py_file':
            module = self.load_module_from_file_path(
                file_path=override['params']['path'])
        elif override['type'] == 'py_module':
            module = self.load_module_from_module_path(
                module_path=override['params']['path'])
        elif override['type'] == 'py_obj':
            module = override['params']['obj']
        elif override['type'] == 'dill_file':
            module = self.load_module_from_dill_path(
                dill_path=override['params']['path'])
        return module

    def load_module_from_file_path(self, file_path=None, module_name=None):
        module_name = module_name or 'random_%s' % random.randint(1, int(1e4))
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def load_module_from_module_path(self, module_path=None):
        return importlib.import_module(module_path)

    def load_module_from_dill_path(self, dill_path=None):
        with open(dill_path, 'rb') as f: return dill.load(f)

    def get_job_module_name(self, job=None, cfg=None):
        module_name = job['job_type']
        return module_name
