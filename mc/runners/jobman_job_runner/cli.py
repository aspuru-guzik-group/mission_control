import argparse
import importlib
import os
import random

import dill

def main():
    parsed_args, command_args = get_cli_args()
    cfg_file_path = parsed_args.config
    cfg = load_cfg_from_file(file_path=cfg_file_path,
                       file_type=parsed_args.config_type)
    job_runner = build_job_runner_per_cfg(cfg=cfg)
    if parsed_args.command == 'tick': job_runner.tick()

def get_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True,
                        type=os.path.expanduser,
                        help=("Path to config file."
                              " File type will be determined by extension,"
                              " or by --config-type argument"))
    parser.add_argument('--config-type')
    parser.add_argument('command')
    return parser.parse_known_args()

def load_cfg_from_file(file_path=None, file_type=None):
    file_type = file_type or os.path.splitext(file_path)[-1]
    if file_type == 'py': cfg = load_module_from_file_path(file_path=file_path)
    elif file_type == 'dill':
        with open(file_path, 'rb') as f: cfg = dill.load(f)
    else: raise Exception("unknown file_type '{file_type}'".format(
        file_type=file_type))
    return cfg

def load_module_from_file_path(file_path=None, module_name=None):
    module_name = module_name or 'random_%s' % random.randint(1, int(1e4))
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def build_job_runner_per_cfg(cfg=None):
    job_runner = getattr(cfg, 'job_runner', None)
    if not job_runner: raise NotImplementedError()
    return job_runner

class JobClient(object):
    def __init__(self, mc_client=None, queue_key=None):
        self.mc_client = mc_client
        self.queue_key = queue_key

    def __getattr__(self, name): return getattr(self.mc_client, name)

    def claim_jobs(self, params=None):
        return self.mc_client.claim_job_queue_items(
            queue_key=self.queue_key, params=params)['items']

if __name__ == '__main__': main()
