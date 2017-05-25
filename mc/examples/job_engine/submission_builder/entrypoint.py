import os
import subprocess
import sys
from unittest.mock import patch

from mc.job_engines.job_engine import JobEngine


def main():
    job_engine = JobEngine()
    import my_job_module
    job = {'job_spec': {'job_type': my_job_module.__name__}}
    submission_meta = job_engine.build_submission(job=job)
    subprocess.run(
        ('cd {dir} && bash {entrypoint}').format(**submission_meta),
        shell=True, check=True
    )
    for log_name in ['stdout', 'stderr']:
        log_path = submission_meta['std_log_files'][log_name]
        with open(log_path) as log: print(log.read())
    print("job_dir contents:\n%s" % "\n".join(
        os.listdir(submission_meta['dir'])))

if __name__ == '__main__':
    this_dir = os.path.dirname(__file__)
    sys_path_w_this_dir = [sys.path[0], this_dir, *sys.path[1:]]
    with patch.object(sys, 'path', new=sys_path_w_this_dir): main()
