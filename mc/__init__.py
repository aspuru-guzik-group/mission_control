VERSION = (0, 0, 1,)
def get_version(): return '.'.join([str(part) for part in VERSION])
__version__ = get_version()
