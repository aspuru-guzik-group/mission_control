import argparse
import os

import yaml
from jobman.jobman import JobMan
from mc.mc_utils.dot_spec_loader import DotSpecLoader
from mc.mc_client import mission_control_client

from . import job_runner

import logging
logging.basicConfig(filename='/tmp/foo.log', level=logging.DEBUG)


class JobClient(object):
    def __init__(self, mc_client=None, queue_key=None):
        self.mc_client = mc_client
        self.queue_key = queue_key

    def __getattr__(self, name): return getattr(self.mc_client, name)

    def claim_jobs(self, params=None):
        return self.mc_client.claim_job_queue_items(
            queue_key=self.queue_key, params=params)['items']

def generate_job_submission_factory(spec=None):
    factory_cls = DotSpecLoader.load_from_dot_spec(spec['type']) 
    return factory_cls.generate_from_spec(spec=spec)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True,
                        type=os.path.expanduser,
                        help='path to config yaml file')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    cfg = yaml.load(open(parsed_args.config))
    if parsed_args.command == 'tick':
        mc_client = mission_control_client.MissionControlClient(
            base_url=cfg['mc_url'])
        job_client = JobClient(mc_client=mc_client,
                               queue_key=cfg['job_queue_key'])
        job_submission_factory = generate_job_submission_factory(
            spec=cfg['job_submission_factory_spec'])
        runner = job_runner.JobRunner(
            job_client=job_client,
            jobman=JobMan(**cfg.get('jobman_cfg', {})),
            job_submission_factory=job_submission_factory,
            logging_cfg=cfg.get('logging_cfg')
        )
        runner.tick()
