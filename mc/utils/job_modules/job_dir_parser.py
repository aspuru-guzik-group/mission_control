from mc.utils.dispatcher import Dispatcher

from . import utils


class JobDirParser(object):
    def __init__(self, dispatcher=None, **kwargs):
        self.dispatcher = dispatcher or Dispatcher()

    @classmethod
    def parse_job_dir(cls, **kwargs):
        return cls(**kwargs)._parse_job_dir(**kwargs)

    def _parse_job_dir(self, job_dir=None):
        job = utils.load_job_dir_component(job_dir=job_dir,
                                           component_name='job_meta')
        return self.dispatcher.dispatch(
            module_name=job['job_type'],
            command='parse_work_dir',
            work_dir=utils.get_job_dir_component_path(
                job_dir=job_dir, component_name='work_dir'),
        )
