import importlib


class ModuleLoader(object):
    """A class that helps with loading job modules."""

    class LoadError(Exception):
        def __init__(self, module_name=None):
            msg = "Could not load module '{module_name}'".format(
                module_name=module_name)
            super().__init__(msg)

    def __init__(self, overrides=None):
        self.overrides = overrides or {}

    @classmethod
    def load_module(cls, module_name=None, overrides=None):
        return cls(overrides=overrides)._load_module(module_name=module_name)

    def _load_module(self, module_name=None):
        """
        Args:
            module_name (str): name of module to load, in the form of a python
                dot path.

        Returns:
            module (module): the loaded module.
        """
        try:
            if module_name in self.overrides:
                return self.overrides[module_name]
            return importlib.import_module(module_name)
        except Exception as exc:
            raise self.LoadError(module_name=module_name) from exc
