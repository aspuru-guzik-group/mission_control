import os
import sys

from mc.flows.flow_engine import FlowEngine

_DIR = os.path.dirname(__file__)

def main():
    sys.path.insert(1, _DIR)
    flow_spec = {
        'label': 'test_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'my_custom_proxying_task',
                'task_params': {
                    'param1': 'value1',
                    'param2': 'value2',
                },
                'precursors': ['ROOT'],
            }
        ]
    }
    flow_engine = FlowEngine()
    flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    try: flow_engine.run_flow(flow=flow)
    except flow_engine.FlowError as exc:
        for error in exc.flow.data.get('errors', []): print(error)
        raise exc

if __name__ == '__main__': main()
