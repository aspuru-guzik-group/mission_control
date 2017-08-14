from mc.utils.dispatcher import Dispatcher

from . import utils


class JobDirParser(object):
    class BuildError(Exception):
        def __init__(self, msg=None, job_dir=None, **kwargs):
            prelude = (
                "Could not parse job_dir; job_dir was '{job_dir}'."
            ).format(job_dir=job_dir)
            msg = (prelude + ' Details: ' + msg) if msg else prelude
            super().__init__(msg, **kwargs)

    def __init__(self, dispatcher=None, **kwargs):
        self.dispatcher = dispatcher or Dispatcher()

    @classmethod
    def parse_job_dir(cls, **kwargs):
        return cls(**kwargs)._parse_job_dir(**kwargs)

    def _parse_job_dir(self, job_dir=None, parse_work_dir_fn=None):
        job = utils.load_job_dir_component(
            job_dir=job_dir, component_name='job_meta')
        try:
            common_parse_kwargs = {
                'job': job,
                'work_dir': utils.get_job_dir_component_path(
                    job_dir=job_dir, component_name='work_dir')
            }
            if parse_work_dir_fn:
                parse_results = parse_work_dir_fn(**common_parse_kwargs)
            else:
                parse_results = self.dispatcher.dispatch(
                    module_name=job['job_type'],
                    command='parse_work_dir',
                    **common_parse_kwargs
                )
            return parse_results
        except Exception as exc:
            raise self.ParseError(
                msg='Failed to parse work_dir.',
                job_dir=self.job_dir
            ) from exc
