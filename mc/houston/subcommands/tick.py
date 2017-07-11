from .base_houston_subcommand import BaseHoustonSubcommand


class TickSubcommand(BaseHoustonSubcommand):
    def _run(self):
        self._tick_flow_runner()
        self._tick_job_runner()
        self._tick_jobman()

    def _tick_flow_runner(self): self.utils.get_flow_runner().tick()
    def _tick_job_runner(self): self.utils.get_job_runner().tick()
    def _tick_jobman(self): self.utils.get_jobman().tick()

Subcommand = TickSubcommand
