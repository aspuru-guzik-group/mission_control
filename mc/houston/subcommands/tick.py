import time
from ._base_houston_subcommand import BaseHoustonSubcommand


class TickSubcommand(BaseHoustonSubcommand):
    DEFAULT_TICK_INTERVAL = 1.0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tick_counter = 0

    def add_arguments(self, parser=None):
        parser.add_argument('--until_finished', action='store_true',
                            default=False)
        parser.add_argument('--interval', type=float, default=1)
        parser.add_argument('--nticks', type=int, default=1)
        parser.add_argument('--max_ticks', type=int)

    def _run(self):
        self.utils.ensure_db()
        self.utils.ensure_queues()
        if self.kwargs.get('until_finished'): self._tick_until_finished()
        else: self._tick_nticks(nticks=self.kwargs['nticks'])

    def _tick_until_finished(self):
        self._tick_while(condition_fn=self.utils.has_unfinished_mc_records)

    def _tick_nticks(self, nticks=None):
        self._tick_while(condition_fn=self._tick_counter_is_lt_nticks)

    def _tick_counter_is_lt_nticks(self):
        return self.tick_counter < self.kwargs['nticks']

    def _tick_while(self, condition_fn=None):
        tick_interval = self._get_tick_interval()
        while condition_fn():
            self.tick_counter += 1
            self._check_max_ticks()
            self._tick()
            self._sleep(tick_interval)

    def _get_tick_interval(self):
        return self.kwargs.get('interval') or self.DEFAULT_TICK_INTERVAL

    def _check_max_ticks(self):
        max_ticks = self.kwargs['max_ticks']
        if max_ticks and self.tick_counter > max_ticks:
            raise Exception("Exceed max_ticks of '%s'" % max_ticks)

    def _tick(self):
        self.logger.info("Tick #{}".format(self.tick_counter))
        self._tick_flow_runner()
        self._tick_job_runner()
        self._tick_jobman()

    def _tick_flow_runner(self): self.utils.flow_runner.tick()
    def _tick_job_runner(self): self.utils.job_runner.tick()
    def _tick_jobman(self): self.utils.jobman.tick()

    def _sleep(self, sleep_time=None): time.sleep(sleep_time)

Subcommand = TickSubcommand
