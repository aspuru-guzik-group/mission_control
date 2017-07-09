import os
import subprocess
import sys
from unittest.mock import patch

from mc.job_module_utils.dispatcher import JobModuleCommandDispatcher


def main():
    dispatcher = JobModuleCommandDispatcher()
    import my_job_module
    job = {'job_type': my_job_module.__name__}
    job_spec = dispatcher.build_jobdir(job=job)
    subprocess.run(
        ('cd {dir} && bash {entrypoint_file_name}').format(**job_spec),
        shell=True, check=True
    )
    for log_name in ['stdout', 'stderr']:
        log_path = os.path.join(job_spec['dir'],
                                job_spec['std_log_file_names'][log_name])
        with open(log_path) as log: print(log.read())
    print("job_dir contents:\n%s" % "\n".join(os.listdir(job_spec['dir'])))

if __name__ == '__main__':
    this_dir = os.path.dirname(__file__)
    sys_path_w_this_dir = [sys.path[0], this_dir, *sys.path[1:]]
    with patch.object(sys, 'path', new=sys_path_w_this_dir): main()
