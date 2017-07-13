import os
import sys
from unittest.mock import patch

from mc.job_module_utils.job_module_command_dispatcher import (
    JobModuleCommandDispatcher)

def main():
    dispatcher = JobModuleCommandDispatcher()
    import my_job_module
    job = {'job_type': my_job_module.__name__}
    job_spec = dispatcher.build_jobdir(job=job)
    dispatcher.run_jobdir(jobdir=job_spec['dir'])

if __name__ == '__main__':
    this_dir = os.path.dirname(__file__)
    sys_path_w_this_dir = [sys.path[0], this_dir, *sys.path[1:]]
    with patch.object(sys, 'path', new=sys_path_w_this_dir): main()
