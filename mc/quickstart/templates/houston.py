#!/usr/bin/env python
import os
import sys


THIS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(THIS_DIR, '.')))
sys.path.append(os.path.abspath(os.path.join(THIS_DIR, '..')))

if __name__ == '__main__':
    def handle_error(err, msg):
        if hasattr(err, 'msg'): err.msg += msg
        else: err.msg = msg
        raise
    
    try:
        from mc.houston.houston import HoustonCommand
    except ImportError as err:
        handle_error(err, "  is 'mc' in your pythonpath?")
    command = HoustonCommand(
        default_cfg_path=os.path.join(THIS_DIR, 'cfgs', 'houston_cfg.py')
    )
    command.run_from_argv(argv=sys.argv)

