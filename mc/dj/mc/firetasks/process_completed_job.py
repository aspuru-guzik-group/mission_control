from fireworks.core.firework import FWAction, FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class ProcessCompletedJobFireTask(FireTaskBase):
    def run_task(self, fw_spec):
        self.process_completed_job_dir(fw_spec['job_dir'])
        return FWAction()

    def process_completed_job_dir(self, job_dir):
        #@TODO: implement!
        print("processing job dir!")
