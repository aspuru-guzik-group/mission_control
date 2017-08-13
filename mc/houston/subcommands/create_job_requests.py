import json

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = ("Create job requests."
            " Will provision new jobs or link existing jobs.")

    def add_arguments(self, parser):
        defaults = self._get_defaults()
        parser.add_argument(
            '--dry_run',
            help=("Don't actually create jobs. Just return # of jobs that"
                  " would be created."),
            action='store_true'
        )

        parser.add_argument(
            '--limit',
            help="Maximum number of requests to create.",
            type=int,
            default=defaults['limit']
        )

        parser.add_argument(
            '--request_dicts',
            help="JSON list of dicts containing request parameters",
            type=json.loads,
            default=defaults['request_dicts']
        )

    def _get_defaults(self):
        return {
            'dry_run': False,
            'limit': 100,
            'request_dicts': [],
        }

    def _run(self):
        self.session.begin_nested()
        requests = self._create_requests()
        if self.parsed_args['dry_run']:
            provisioning_tallies = {'DRY_RUN': '<nothing created>'}
        else:
            provisioning_tallies = self._provision_jobs_for_requests(requests)
        if self.parsed_args['dry_run']:
            self.session.rollback()
        else:
            self.session.commit()
        return {
            **provisioning_tallies,
            'num_requests_created': len(requests),
        }

    def _create_requests(self):
        requests = [
            self.Request(request_type='job', **request_dict)
            for request_dict in self.parsed_args['request_dicts']
        ]
        return requests

    @property
    def Request(self): return self.db.models.Request

    def _provision_jobs_for_requests(self, requests=None):
        tallies = {'num_existing_jobs_linked': 0, 'num_jobs_created': 0}
        if not requests:
            return tallies
        jobs_by_request = {
            request: self.Job(
                module=request.params['job_module'],
                params=request.params['job_params']
            )
            for request in requests
        }
        jobs_by_hash = {job.hash: job for job in jobs_by_request.values()}
        existing_jobs = (
            self.session.query(self.Job)
            .filter(self.Job.hash.in_(jobs_by_hash.keys()))
            .all()
        )
        existing_jobs_by_hash = {job.hash: job for job in existing_jobs}
        for request, job in jobs_by_request.items():
            if job.hash in existing_jobs_by_hash:
                job = existing_jobs_by_hash[job.hash]
                tallies['num_existing_jobs_linked'] += 1
            else:
                job = jobs_by_hash[job.hash]
                tallies['num_jobs_created'] += 1
            request.child_nodes.append(job)
            request.status = 'PROVISIONED'
        self.session.add_all(requests)
        return tallies

    @property
    def Job(self): return self.db.models.Job
