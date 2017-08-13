import glob
import itertools
import traceback

from a2g2_v2.utils.job_modules import utils as _jmu_utils
from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Build dirs for pending jobs."

    def add_arguments(self, parser):
        defaults = self._get_defaults()
        parser.add_argument(
            '--search_globs',
            help=("Specify 1 or more globs to use for searching for job dirs."
                  " Separate globs with spaces."),
            dest='search_globs',
            nargs='*',
            default=defaults['search_globs']
        )

        parser.add_argument(
            '--dry_run',
            help=("Don't actually process executed job dirs. Just return dirs"
                  " that would have been processed."),
            action='store_true'
        )

        parser.add_argument(
            '--limit',
            help="Maximum number of dirs to process.",
            type=int,
            default=defaults['limit']
        )

    def _get_defaults(self):
        return {
            'dry_run': False,
            'limit': 100,
            'search_globs': [
                (self.houston.utils.job_dirs['executed'] + '/*')
            ],
        }

    def _run(self):
        self.session.begin_nested()
        process_results = self._process_paths(
            paths=self._get_paths_to_process())
        if self.parsed_args['dry_run']:
            self.session.rollback()
        else:
            self.session.commit()
        return process_results

    def _get_paths_to_process(self):
        return itertools.chain(*[
            glob.iglob(search_glob)
            for search_glob in self.parsed_args['search_globs']
        ])

    def _process_paths(self, paths=None):
        results = {'visited': set(), 'unrecognized': [], 'completed': [],
                   'processed': [], 'failed': {}, 'failed_to_process': {}}
        limit = self.parsed_args['limit']
        for path in paths:
            if path in results['visited']:
                continue
            results['visited'].add(path)
            try:
                if _jmu_utils.is_job_dir(path):
                    evaluation = self._process_job_dir(job_dir=path)
                    results['processed'].append(path)
                    if evaluation['status'] == 'COMPLETED':
                        results['completed'].append(path)
                    elif evaluation['status'] == 'FAILED':
                        results['failed'][path] = evaluation['error']
                else:
                    results['unrecognized'].append(path)
            except Exception as exc:
                results['failed_to_process'][path] = traceback.format_exc()
            if limit and len(results['processed']) > limit:
                break
        return results

    def _process_job_dir(self, job_dir=None):
        job = self._get_job_for_job_dir(job_dir)
        evaluation = _jmu_utils.evaluate_completed_job_dir(job_dir)
        job.status = evaluation['status']
        if job.status == 'FAILED':
            job.error = evaluation.get('error', '<unknown error>')
        elif job.status == 'COMPLETED':
            if not self.parsed_args['dry_run']:
                job.archive_meta = self.houston.utils.archiver.ingest(
                    src=job_dir)
        self._update_requests_related_to_job(job)
        self.session.add(job)
        return evaluation

    def _get_job_for_job_dir(self, job_dir=None):
        job_key = _jmu_utils.load_job_dir_component(job_dir=job_dir,
                                                    component_name='job_key')
        return (
            self.session.query(self.db.models.Job)
            .filter_by(key=job_key.strip())
            .first()
        )

    def _update_requests_related_to_job(self, job=None):
        parent_requests = [
            request for request in job.parent_nodes_of_type('Request')
        ]
        for parent_request in parent_requests:
            parent_request.status = job.status
        self.session.add_all(parent_requests)
