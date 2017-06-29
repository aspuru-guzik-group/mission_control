import os

from mc.flows.flow_engine import FlowEngine

_DIR = os.path.dirname(__file__)

def main():
    flow_spec = {
        'label': 'test_flow',
        'tasks': [
            {
                'key': 'spreader',
                'task_type': 'spread',
                'task_params': {
                    'container_task_type': 'mc.tasks.inline_flow',
                    'items': ['item_%s' % i for i in range(3)],
                    'mapping_params': {
                        'skeleton': {
                            'task_type': 'print',
                            'task_params': {'msg': 'TO BE WIRED'}
                        },
                        'wirings': [
                            'ctx.item.value => ctx.skeleton.task_params.msg'
                        ]
                    }
                }
            }
        ]
    }
    flow_engine = FlowEngine()
    flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    flow_engine.run_flow(flow=flow)

if __name__ == '__main__': main()
