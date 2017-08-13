import json
from pathlib import Path
import traceback

from mc.utils.job_modules import constants as _jmu_constants
from mc.utils import dot_spec_loader
from mc.utils import hash_utils

from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    help = "Create new jobs or link existing jobs for pending job requests."

    class InvalidArgsError(Exception):
        pass

    class ParseError(Exception):
        pass

    def add_arguments(self, parser):
        defaults = self._get_defaults()
        parser.add_argument(
            '--request_tag',
            help=("Only parse requests that have this tag.")
        )

        parser.add_argument(
            '--dry_run',
            help=("If true don't write changes to DB."),
            action='store_true'
        )

        parser.add_argument(
            '--limit',
            help="Maximum number of requests to parse.",
            type=int,
            default=defaults['limit']
        )

        parser.add_argument(
            '--parser',
            help="dotspec for fn that parses workdir",
            required=True,
        )

        parser.add_argument(
            '--parser_params',
            help="JSON dict of parser params",
            type=json.loads,
            default=defaults['parser_params'],
        )

        parser.add_argument(
            '--tags_to_apply',
            help=("apply this tag to results. Can specify multiple tags"
                  " separated by spaces"),
            nargs='*',
            dest='tags_to_apply',
            default=defaults['tags_to_apply']
        )

    def _get_defaults(self):
        return {
            'dry_run': False,
            'limit': 100,
            'request_tag': None,
            'tags_to_apply': [],
            'parser_params': {},
        }

    def _run(self):
        self.session.begin_nested()
        parse_results = self._parse_job_requests(
            job_requests=self._get_unparsed_job_requests(),
            parse_fn=self._get_parse_fn()
        )
        if self.parsed_args['dry_run']:
            self.session.rollback()
        else:
            self.session.commit()
        return parse_results

    def _get_parse_fn(self):
        try:
            return dot_spec_loader.load_from_dot_spec(
                dot_spec=self.parsed_args['parser'])
        except Exception as exc:
            raise self.InvalidArgsError("Could not load parser") from exc

    def _get_unparsed_job_requests(self):
        unparsed_job_request_items = (
            self.houston.utils.request_selector.get_items(
                having_request_tag=self.parsed_args['request_tag'],
                having_status='COMPLETED',
                sans_downstream_requests_w_tags=[self._get_parse_request_tag()]
            )
        )
        return [item['value'] for item in unparsed_job_request_items]

    def _get_parse_request_tag(self):
        if not hasattr(self, '_parse_request_tag'):
            parsing_hash = hash_utils.hash_obj(
                [self.parsed_args['parser'], self.parsed_args['parser_params']]
            )
            self._parse_request_tag = 'PARSE:' + parsing_hash
        return self._parse_request_tag

    def _parse_job_requests(self, job_requests=None, parse_fn=None):
        results = {'completed': set(), 'failed': {}}
        for job_request in job_requests:
            try:
                self._parse_job_request(job_request=job_request,
                                        parse_fn=parse_fn)
                results['completed'].add(job_request.key)
            except Exception:
                error = traceback.format_exc()
                results['failed'][job_request.key] = error
                msg = (
                    "failed to parse job_request '{key}', error: {error}"
                ).format(key=job_request.key, error=error)
                self.logger.warn(msg)
        return results

    def _parse_job_request(self, job_request=None, parse_fn=None):
        with self.session.begin_nested():
            parse_request = self._create_parse_request(job_request=job_request)
            self.session.add(parse_request)
        try:
            work_dir = self._get_work_dir_for_job_request(job_request)
            parse_result = parse_fn(
                work_dir=work_dir, **self.parsed_args['parser_params'])
            result = self._process_actions(
                actions=parse_result.get('actions', []),
                decorations=self._get_decorations_for_job_request(job_request)
            )
            parse_request.status = 'COMPLETED'
            return result
        except Exception as exc:
            parse_request.status = 'FAILED'
            parse_request.error = traceback.format_exc()
            self.session.add(parse_request)
            raise exc

    def _create_parse_request(self, job_request=None):
        return self.Request(
            request_type='parse',
            request_tag=self._get_parse_request_tag(),
            instance_key=job_request.key,
            status='PARSING'
        )

    @property
    def Request(self): return self.db.models.Request

    def _get_work_dir_for_job_request(self, job_request=None):
        job = job_request.child_nodes_of_type('Job')[0]
        materialized_job_dir = self.houston.utils.archiver.materialize_as_path(
            meta=job.artifact_meta)
        return Path(
            materialized_job_dir,
            _jmu_constants.JOB_DIR_COMPONENT_PATHS['work_dir']
        )

    def _get_decorations_for_job_request(self, job_request=None):
        decorations = {}
        tags_to_apply = self.parsed_args['tags_to_apply']
        if tags_to_apply:
            decorations['tags'] = tags_to_apply
        ancestor_keys = job_request.params.get('ancestor_keys')
        if ancestor_keys:
            decorations['ancestor_keys'] = ancestor_keys
        return decorations

    def _process_actions(self, actions=None, decorations=None):
        results = []
        with self.session.begin_nested():
            for action in actions:
                result = self._process_action(action, decorations=decorations)
                results.append(result)
        return results

    def _process_action(self, action=None, decorations=None):
        if action['type'] == 'upsert':
            return self._process_upsert_action(action=action,
                                               decorations=decorations)

    def _process_upsert_action(self, action=None, decorations=None):
        extra_actions = []
        tags = decorations.get('tags')
        if tags:
            extra_actions.append(('tags', '$addToSet', tags))
        ancestor_keys = decorations.get('ancestor_keys')
        if ancestor_keys:
            extra_actions.append(('add_ancestors_by_key', '$call', (),
                                  {'keys': ancestor_keys}))
        action['params']['updates'].extend(extra_actions)
        return self.db.execute_action(action=action)
