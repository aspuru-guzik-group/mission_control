import json

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Get an existing job, or create a new job if it does not yet exist"

    def add_arguments(self, parser=None):
        parser.add_argument('--job_type', required=True)
        parser.add_argument('--job_params', type=json.loads, required=True)

    def _run(self):
        job = self.Job(
            job_type=self.parsed_args['job_type'],
            job_params=self.parsed_args['job_params'],
        )
        existing_job = self._get_existing_job(job=job)
        if existing_job:
            job = existing_job
        if not existing_job:
            job = self._create_job(job=job)
        return job.to_dict()

    @property
    def Job(self): return self.houston.utils.db.models.Job

    def _get_existing_job(self, job=None):
        return (
            self.session.query(self.Job)
            .filter_by(job_hash=job.job_hash)
            .first()
        )

    @property
    def session(self): return self.houston.utils.db.session

    def _create_job(self, job=None):
        with self.session.begin(subtransactions=True):
            self.db.session.add(job)
        return job
