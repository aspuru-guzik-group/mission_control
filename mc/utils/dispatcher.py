import logging

from .module_loader import ModuleLoader


class Dispatcher(object):
    def __init__(self, load_module_fn=None, logger=None):
        """
        Args:
            load_module_fn (fn): A fn to use for loading modules. See
                :meth:`_default_load_module_fn for more info.
            logger: (logging.Logger): a logger.
        """
        self.logger = logger or logging
        self.load_module_fn = load_module_fn or self._default_load_module_fn

    def _default_load_module_fn(self, module_name=None):
        """Load a job module.

        Args:
            module_name (str): name of the module to load

        Returns:
            module (module): loaded module
        """
        return ModuleLoader().load_module(module_name=module_name)

    def dispatch(self, module_name=None, command=None, **kwargs):
        module = self.load_module_fn(module_name=module_name)
        fn_for_command = getattr(module, command)
        return fn_for_command(**kwargs)
