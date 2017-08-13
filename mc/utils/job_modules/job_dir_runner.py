import traceback

from mc.utils.dispatcher import Dispatcher

from . import utils


class JobDirRunner(object):
    def __init__(self, dispatcher=None):
        self.dispatcher = dispatcher or Dispatcher()

    def run_job_dir(self, job_dir=None, fake=None):
        job = utils.load_job_dir_component(job_dir=job_dir,
                                           component_name='job_meta')
        try:
            self.dispatcher.dispatch(
                module_name=job['module'],
                command='run_work_dir',
                work_dir=utils.get_job_dir_component_path(
                    job_dir=job_dir, component_name='work_dir'
                ),
                fake=fake
            )
            completed_checkpoint_path = utils.get_job_dir_component_path(
                job_dir=job_dir, component_name='completed_checkpoint')
            completed_checkpoint_path.touch()
        except Exception as exc:
            failure_checkpoint_path = utils.get_job_dir_component_path(
                job_dir=job_dir, component_name='failure_checkpoint')
            with failure_checkpoint_path.open(mode='w') as f:
                f.write(traceback.format_exc())
