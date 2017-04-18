from mc.task_handlers.base_task_handler import BaseTaskHandler

from ..parse.qchem_computation_parser import QChemComputationParser

class TaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None, **kwargs):
        submission = task_context['flow_ctx']['submission']
        input_dir = submission['inputs_dir'] + '/input_dir'
        parser = QChemComputationParser(input_dir=input_dir)
        computation_meta = parser.extract_computation_meta(job_dir=input_dir)
        task['data']['computation_meta'] = computation_meta
        task['status'] = 'COMPLETED'
