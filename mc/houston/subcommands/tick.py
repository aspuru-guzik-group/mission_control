import time
from ._base_houston_subcommand import BaseHoustonSubcommand


class TickSubcommand(BaseHoustonSubcommand):
    DEFAULT_TICK_INTERVAL = 1.0

    def add_arguments(self, parser=None):
        parser.add_argument('--until_finished', action='store_true',
                            default=False)
        parser.add_argument('--tick_interval', type=float, default=1)

    def _run(self):
        self.tick_counter = 0
        self.utils.ensure_db()
        self.utils.ensure_queues()
        if self.kwargs.get('until_finished'): self._tick_until_finished()
        else: self._tick()

    def _tick_until_finished(self):
        tick_interval = self._get_tick_interval()
        while self.utils.has_unfinished_mc_records():
            self._tick()
            time.sleep(tick_interval)
        self.logger.info('Finished!')

    def _get_tick_interval(self):
        return self.kwargs.get('tick_interval') or \
                self.cfg.get('TICK_INTERVAL', None) or \
                self.DEFAULT_TICK_INTERVAL

    def _tick(self):
        self.tick_counter += 1
        self.logger.info("Tick #{}".format(self.tick_counter))
        self._tick_flow_runner()
        self._tick_job_runner()
        self._tick_jobman()

    def _tick_flow_runner(self): self.utils.flow_runner.tick()
    def _tick_job_runner(self): self.utils.job_runner.tick()
    def _tick_jobman(self): self.utils.jobman.tick()

Subcommand = TickSubcommand
