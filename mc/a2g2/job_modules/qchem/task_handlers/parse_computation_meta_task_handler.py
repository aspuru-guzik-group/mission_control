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
            computation_meta=computation_meta,
            task=task,
            task_context=task_context)
        task['status'] = 'COMPLETED'

    def generate_chemthing_actions(self, computation_meta=None, task=None,
                                   task_context=None):
        chemthing_actions = [
            {
                'key': computation_meta['uuid'],
                'updates': {
                    'props': {
                        'a2g2:prop:artifacts': (
                            self.get_artifacts_from_task_context(
                                task_context=task_context)),
                        'a2g2:prop:command_meta': (
                            computation_meta['command_meta']),
                        'a2g2:prop:execution_meta': (
                            computation_meta['execution_meta']),
                    },
                    'precursors': self.get_precursors(task=task)
                }
            }
        ]
        return chemthing_actions

    def get_artifacts_from_task_context(self, task_context=None):
        return task_context['flow_ctx']['job']['job_spec'].get(
            'inputs', {}).get('artifacts')

    def get_precursors(self, task=None):
        precursors = task['task_params'].get('precursors', [])
        precursor_uuid_dict = {}
        for precursor in precursors:
            precursor_uuid = precursor.get('uuid')
            if precursor_uuid: precursor_uuid_dict[precursor_uuid] = True
        return precursor_uuid_dict
