from mc.flows.flow_engine import FlowEngine

def main():
    flow_spec = {
        'label': 'example_flow',
        'tasks': [
            {
                'key': 'task_1',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_1.'},
            },
            {
                'key': 'task_2',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_2.'},
            },
            {
                'key': 'task_3',
                'task_type': 'print',
                'task_params': {'msg': 'I am task_3.'},
            }
        ]
    }
    flow_engine = FlowEngine()
    flow = flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
    flow_engine.run_flow(flow=flow)

if __name__ == '__main__': main()
