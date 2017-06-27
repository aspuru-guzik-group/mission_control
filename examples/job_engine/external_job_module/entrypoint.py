import os
import sys
from unittest.mock import patch

from mc.job_engines.job_engine import JobEngine

def main():
    job_engine = JobEngine()
    import my_job_module
    job = {'job_spec': {'job_type': my_job_module.__name__}}
    submission_meta = job_engine.build_job_submission(job=job)
    job_engine.run_job_submission(submission_dir=submission_meta['dir'])

if __name__ == '__main__':
    this_dir = os.path.dirname(__file__)
    sys_path_w_this_dir = [sys.path[0], this_dir, *sys.path[1:]]
    with patch.object(sys, 'path', new=sys_path_w_this_dir): main()
