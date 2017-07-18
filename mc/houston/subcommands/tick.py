import argparse
import json
import os
import time

from ._base_houston_subcommand import BaseHoustonSubcommand


class TickSubcommand(BaseHoustonSubcommand):
    DEFAULT_TICK_INTERVAL = 1.0
    TICKEE_NAMES = ['flow', 'flow_runner', 'job_runner', 'jobman', 'all']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tick_counter = 0

    def add_arguments(self, parser=None):
        parser.add_argument(
            'tickee', help="The component to tick",
            choices=self.TICKEE_NAMES
        )
        parser.add_argument('--interval', type=float, default=1)
        parser.add_argument('--nticks', type=int, default=1)
        parser.add_argument('--max_ticks', type=int)

    def _run(self):
        output_for_tickee = self._dispatch_on_tickee()
        print(output_for_tickee)

    def _dispatch_on_tickee(self):
        tickee_name = self.parsed_args['tickee']
        fn_for_tickee = self._get_fn_for_tickee(tickee_name=tickee_name)
        output_for_tickee = fn_for_tickee()
        return output_for_tickee

    def _get_fn_for_tickee(self, tickee_name):
        fn_for_tickee = getattr(self, '_run_for_%s_tickee' % tickee_name)
        return fn_for_tickee

    def _run_for_flow_tickee(self):
        self._update_args_for_flow_tickee()
        return self._run_for_flow_spec()

    def _update_args_for_flow_tickee(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--indent', type=int)
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--flow_spec', type=json.loads,
                           help="A flow_spec JSON string")
        group.add_argument('--flow_spec_file', help="A flow spec JSON file",
                           type=self._json_path)
        parser.add_argument('--until_finished', action='store_true')
        parsed, unparsed = parser.parse_known_args(self.unparsed_args)
        self._update_args(parsed=parsed, unparsed=unparsed)

    def _json_path(self, path):
        with open(os.path.expanduser(path)) as f: return json.load(f)

    def _update_args(self, parsed=None, unparsed=None):
        self.parsed_args.update(vars(parsed))
        self.unparsed_args = unparsed

    def _run_for_flow_spec(self):
        flow_spec = self.parsed_args.get('flow_spec') or \
                self.parsed_args.get('flow_spec_file')
        flow = self.utils.flow_engine.flow_spec_to_flow(flow_spec=flow_spec)
        output_for_flow = self._run_for_flow(flow=flow)
        return output_for_flow

    def _run_for_flow(self, flow=None):
        flow_ref = {'flow': flow}
        def _get_flow():
            nonlocal flow_ref
            return flow_ref['flow']

        def _set_flow(updated_flow):
            nonlocal flow_ref
            flow_ref['flow'] = updated_flow

        self._tick_while(
            tick_fn=self._get_tick_fn_for_flow(
                get_flow=_get_flow, set_flow=_set_flow),
            condition_fns=self._get_condition_fns_for_flow(get_flow=_get_flow)
        )
        output_for_flow = self._get_output_for_flow(flow=_get_flow())
        return output_for_flow

    def _get_tick_fn_for_flow(self, get_flow=None, set_flow=None):
        def _tick_fn():
            flow = get_flow()
            self.utils.flow_engine.tick_flow(flow=get_flow())
            set_flow(flow)
        return _tick_fn

    def _get_condition_fns_for_flow(self, get_flow=None):
        condition_fns = []
        def _flow_is_unfinished():
            return get_flow().status not in {'COMPLETED', 'FAILED'}
        condition_fns.append(_flow_is_unfinished)
        if not self.parsed_args.get('until_finished') and \
           self.parsed_args.get('nticks'):
            condition_fns.append(self._tick_counter_is_lt_nticks)
        return condition_fns

    def _tick_while(self, tick_fn=None, condition_fns=None):
        condition_fns = condition_fns or []
        tick_interval = self._get_tick_interval()
        while all([fn() for fn in condition_fns]):
            self.tick_counter += 1
            if self.parsed_args.get('verbose'):
                self.logger.info("{self} Tick #{tick_counter}".format(
                    self=self, tick_counter=self.tick_counter))
            self._check_max_ticks()
            tick_fn()
            self._sleep(tick_interval)

    def _get_common_condition_fns(self):
        common_condition_fns = []
        if self.parsed_args.get('until_finished'):
            common_condition_fns.append(self._has_unfinished_mc_records)
        elif self.parsed_args.get('nticks'):
            common_condition_fns.append(self._tick_counter_is_lt_nticks)
        return common_condition_fns

    def _has_unfinished_mc_records(self):
        return self.utils.has_unfinished_mc_records()

    def _tick_counter_is_lt_nticks(self):
        return self.tick_counter < self.parsed_args['nticks']

    def _get_tick_interval(self):
        return self.parsed_args.get('interval') or self.DEFAULT_TICK_INTERVAL

    def _check_max_ticks(self):
        max_ticks = self.parsed_args['max_ticks']
        if max_ticks and self.tick_counter > max_ticks:
            raise Exception("Exceed max_ticks of '%s'" % max_ticks)

    def _get_output_for_flow(self, flow=None):
         return json.dumps(self.utils.flow_engine.flow_to_flow_dict(flow=flow),
                           indent=self.parsed_args.get('indent'))

    def _run_for_flow_runner_tickee(self):
        self._ensure_mc()
        self._update_args_for_flow_runner_tickee()
        self._tick_while(
            tick_fn=self._tick_flow_runner,
            condition_fns=self._get_common_condition_fns()
        )

    def _ensure_mc(self):
        self.utils.ensure_db()
        self.utils.ensure_queues()

    def _update_args_for_flow_runner_tickee(self):
        parser = argparse.ArgumentParser()
        self._add_mc_runner_arguments(parser=parser)
        parsed, unparsed = parser.parse_known_args(self.unparsed_args)
        self._update_args(parsed=parsed, unparsed=unparsed)

    def _add_mc_runner_arguments(self, parser=None):
        parser.add_argument('--until_finished', action='store_true')

    def _run_for_job_runner_tickee(self):
        self._ensure_mc()
        self._update_args_for_job_runner_tickee()
        self._tick_while(
            tick_fn=self._tick_job_runner,
            condition_fns=self._get_common_condition_fns()
        )

    def _update_args_for_job_runner_tickee(self):
        parser = argparse.ArgumentParser()
        self._add_mc_runner_arguments(parser=parser)
        parsed, unparsed = parser.parse_known_args(self.unparsed_args)
        self._update_args(parsed=parsed, unparsed=unparsed)

    def _run_for_jobman_tickee(self):
        self._update_args_for_jobman_tickee()
        self._tick_while(
            tick_fn=self._tick_jobman,
            condition_fns=self._get_common_condition_fns()
        )

    def _update_args_for_jobman_tickee(self): pass

    def _run_for_all_tickee(self):
        for update_args_fn in [
            self._update_args_for_flow_runner_tickee,
            self._update_args_for_job_runner_tickee,
            self._update_args_for_jobman_tickee,
        ]: update_args_fn()
        self._tick_while(
            tick_fn=self._tick_all,
            condition_fns=self._get_common_condition_fns()
        )

    def _tick_all(self):
        self._tick_flow_runner()
        self._tick_job_runner()
        self._tick_jobman()

    def _tick_flow_runner(self): self.utils.flow_runner.tick()
    def _tick_job_runner(self): self.utils.job_runner.tick()
    def _tick_jobman(self): self.utils.jobman.tick()

    def _sleep(self, sleep_time=None): time.sleep(sleep_time)

Subcommand = TickSubcommand
