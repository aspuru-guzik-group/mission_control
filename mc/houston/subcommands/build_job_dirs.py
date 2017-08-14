from pathlib import Path
import time
import traceback

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Build dirs for pending jobs."

    def add_arguments(self, parser):
        defaults = self._get_defaults()
        parser.add_argument(
            '--tag',
            help=("Only build jobs that have this tag. Can specify"
                  " this arg multiple times to filter for multiple tags."),
            dest='tags',
            action='append'
        )

        parser.add_argument(
            '--parent_request_tag',
            help=("Only build jobs which have a parent request with this tag."
                  " Can specify this arg multiple times to filter for multiple"
                  " tags."),
            dest='parent_request_tags',
            action='append'
        )

        parser.add_argument(
            '--output_dir',
            help=("Where to put created job dirs. If dir does not exist it"
                  " will be created."),
            default=defaults['output_dir']
        )

        parser.add_argument(
            '--dry_run',
            help=("Don't actually build jobs. Just return # of jobs that would"
                  " be built."),
            action='store_true'
        )

        parser.add_argument(
            '--limit',
            help="Maximum number of jobs to claim.",
            type=int,
            default=defaults['limit']
        )

        parser.add_argument(
            '--job_dir_tpl',
            help=("Template for naming job_dirs"),
            default=defaults['job_dir_tpl']
        )

    def _get_defaults(self):
        return {
            'dry_run': False,
            'limit': 100,
            'job_dir_tpl': '{timestamp}.{key}',
            'output_dir': self.houston.utils.job_dirs['pending'],
            'parent_request_tags': [],
            'tags': [],
        }

    def _run(self):
        self.session.begin_nested()
        output_dir = self._get_output_dir()
        claimed_jobs = self._get_and_claim_jobs()
        if self.parsed_args['dry_run']:
            build_results = {'DRY_RUN': '<nothing built>'}
        else:
            build_results = self._build_jobs(jobs=claimed_jobs,
                                             output_dir=output_dir)
        if self.parsed_args['dry_run']:
            self.session.rollback()
        else:
            self.session.commit()
        if build_results.get('errors'):
            raise Exception("\n".join(build_results['errors']))
        return {
            **build_results,
            'output_dir': output_dir,
            'num_claimed': len(claimed_jobs),
        }

    def _get_output_dir(self):
        output_dir = self.parsed_args['output_dir']
        if not self.parsed_args['dry_run']:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        return output_dir

    def _get_and_claim_jobs(self):
        claim_time = time.time()
        modified_clause = (self.Job.modified < claim_time)
        q = (
            self.session.query(self.Job)
            .filter(self.Job.status == 'PENDING')
            .filter(modified_clause)
        )
        for tag_name in self.parsed_args['tags']:
            q = q.filter(self.Job.tags_set.any(name=tag_name))
        if self.parsed_args['limit']:
            q = q.limit(self.parsed_args['limit'])
        if self.parsed_args['parent_request_tags']:
            q = (
                q.join(self.Job.parents, aliased=True)
                .join(self.Request, aliased=True, from_joinpoint=True)
                .filter(
                    self.Request.request_tag.in_(
                        self.parsed_args['parent_request_tags']
                    )
                )
                .reset_joinpoint()
            )
        jobs_to_claim = q.all()
        if not jobs_to_claim:
            return []
        keys_clause = self.Job.key.in_([j.key for j in jobs_to_claim])
        (
            self.session.query(self.Job)
            .filter(keys_clause)
            .filter(modified_clause)
            .update(
                {'status': 'BUILDING', 'modified': claim_time},
                synchronize_session='fetch'
            )
        )
        claimed_jobs = (
            self.session.query(self.Job)
            .filter(keys_clause)
            .filter(self.Job.modified == claim_time)
            .all()
        )
        return claimed_jobs

    @property
    def Job(self): return self.db.models.Job

    @property
    def Request(self): return self.db.models.Request

    def _build_jobs(self, jobs=None, output_dir=None):
        num_built = 0
        errors = []
        job_dir_builder = self._get_job_dir_builder()
        for job in jobs:
            try:
                self._build_job(job=job,
                                job_dir_builder=job_dir_builder,
                                output_dir=output_dir)
                num_built += 1
            except:
                errors.append(traceback.format_exc())
        with self.session.begin(subtransactions=True):
            self.session.add_all(jobs)
        return {'num_built': num_built, 'errors': errors}

    def _get_job_dir_builder(self):
        from mc.utils.job_modules.job_dir_builder import JobDirBuilder
        return JobDirBuilder()

    def _build_job(self, job=None, job_dir_builder=None,
                   output_dir=None):
        try:
            job_dir_builder.build_job_dir(
                job_dict=job.to_dict(),
                output_dir=self._get_job_output_dir(job=job,
                                                    parent_dir=output_dir)
            )
            job.status = 'BUILT'
        except Exception as exc:
            job.status = 'FAILED'
            raise

    def _get_job_output_dir(self, job=None, parent_dir=None):
        tpl = self.parsed_args.get('job_dir_tpl')
        job_dir_name = tpl.format(timestamp=int(time.time()), key=job.key)
        return Path(parent_dir, job_dir_name)
