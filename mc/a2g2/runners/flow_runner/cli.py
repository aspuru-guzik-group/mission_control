import argparse

from mc.mc_client import mission_control_client
from . import flow_runner

import logging
logging.basicConfig(filename='/tmp/flow_runner.log', level=logging.DEBUG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mc-url', help='url of MissionControl server')
    parser.add_argument('command')
    parsed_args, command_args = parser.parse_known_args()
    if parsed_args.command == 'tick':
        mc_client = mission_control_client.MissionControlClient(
            base_url=parsed_args.mc_url)
        runner = flow_runner.FlowRunner(mc_client=mc_client)
        runner.tick()
