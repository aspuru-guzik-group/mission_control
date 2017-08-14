import importlib
import os
import random


def load_module_from_path(path=None, module_name=None):
    module_name = module_name or get_random_module_name_for_path(path)
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def get_random_module_name_for_path(path=None):
    return '{random_prefix}_{basename}'.format(
        random_prefix=random.randint(0, 1e6),
        basename=os.path.basename(os.path.splitext(path)[0])
    )
