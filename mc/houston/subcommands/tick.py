from .base_houston_subcommand import BaseHoustonSubcommand


class TickSubcommand(BaseHoustonSubcommand):
    def _run(self):
        self.utils.ensure_db()
        self.utils.ensure_queues()
        self._tick_flow_runner()
        self._tick_job_runner()
        self._tick_jobman()

    def _tick_flow_runner(self): self.utils.flow_runner.tick()
    def _tick_job_runner(self): self.utils.job_runner.tick()
    def _tick_jobman(self): self.utils.jobman.tick()

Subcommand = TickSubcommand
