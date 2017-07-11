from .base_houston_subcommand import BaseHoustonSubcommand

class FlushMcRecordsSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        parser.add_argument('--record_type', required=True)

    def _run(self):
        self.utils.ensure_db()
        record_type = self.kwargs['record_type']
        if record_type == 'ALL': self._flush_all()
        elif record_type == 'Flow': self._flush_flows()
        elif record_type == 'Job': self._flush_jobs()
        elif record_type == 'Queue': self._flush_queues()
        else:
            raise Exception("Unknown record type '{record_type}'".format(
                record_type=record_type))

    def _flush_all(self): self.utils.mc_dao.flush_mc_db()
    def _flush_flows(self): self.utils.mc_dao.delete_items(item_type='Flow')
    def _flush_jobs(self): self.utils.mc_dao.delete_items(item_type='Job')
    def _flush_queues(self): self.utils.mc_dao.delete_items(item_type='Queue')

Subcommand = FlushMcRecordsSubcommand
