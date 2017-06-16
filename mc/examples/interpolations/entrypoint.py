import os

from mc.flow_engines.flow_engine import FlowEngine

_DIR = os.path.dirname(__file__)

def main():
    flow_spec = {
        'label': 'test_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'wire',
                'task_params': {
                    'wirings': [
                        {'dest': 'ctx.flow.tasks.task_2.task_params.msg',
                         'value': 'msg set by task_1'}
                    ]
                },
                'data': {'some_key': 'value from task_1.data'},
            },
            {
                'key': 'task_2',
                'task_type': 'print',
                'task_params': {'msg': 'I will be set from task_1'},
            },
            {
                'key': 'task_3',
                'task_type': 'print',
                'task_params': {'msg': '$ctx.flow.tasks.task_1.data.some_key'},
            }
        ]
    }
    flow_engine = FlowEngine()
    flow = flow_engine.generate_flow(flow_spec=flow_spec)
    flow_engine.run_flow(flow=flow)

if __name__ == '__main__': main()
