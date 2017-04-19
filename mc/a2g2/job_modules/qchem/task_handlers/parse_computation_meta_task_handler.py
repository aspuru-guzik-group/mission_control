from mc.task_handlers.base_task_handler import BaseTaskHandler

from ..parse.qchem_computation_parser import QChemComputationParser

class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        submission = task_context['flow_ctx']['submission']
        input_dir = submission['inputs_dir'] + '/input_dir'
        parser = QChemComputationParser(input_dir=input_dir)
        computation_meta = parser.extract_computation_meta(job_dir=input_dir)
        task['data']['computation_meta'] = computation_meta
        task['data']['chemthing_actions'] = self.generate_chemthing_actions(
            computation_meta=computation_meta)
        task['status'] = 'COMPLETED'

    def generate_chemthing_actions(self, computation_meta=None):
        chemthing_actions = [
            {
                'key': computation_meta['uuid'],
                'updates': {
                    'props': {
                        'a2g2:prop:artifacts': 'HERE!!!',
                        'a2g2:prop:command_meta': (
                            computation_meta['command_meta']),
                        'a2g2:prop:execution_meta': (
                            computation_meta['execution_meta']),
                    }
                }
            }
        ]
        return chemthing_actions
