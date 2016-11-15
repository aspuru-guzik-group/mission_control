from fireworks.core.firework import Firework, FWAction, FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from .process_odyssey_job_request import ProcessOdysseyJobRequestFireTask

@explicit_serialize
class ProcessJobRequestFireTask(FireTaskBase):
   def run_task(self, fw_spec):
       print("ProcessJobRequestFireTask")
       # Determines queue to allocate job request to.
       # @TODO: implement! For now, just assigns to odyssey queue.
       ftask = ProcessOdysseyJobRequestFireTask()
       fwork = Firework([ftask])
       fwork.spec['_category'] = 'odyssey_job_requests'
       return FWAction(additions=[fwork])
