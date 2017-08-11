from ._base_houston_subcommand import BaseHoustonSubcommand


class FlushMcRecordsSubcommand(BaseHoustonSubcommand):
    def add_arguments(self, parser=None):
        parser.add_argument('--record_type', required=True)

    def _run(self):
        self.utils.ensure_db()
        record_type = self.parsed_args['record_type']
        if record_type == 'ALL':
            self._flush_all()
        elif record_type == 'flow':
            self._flush_flows()
        elif record_type == 'job':
            self._flush_jobs()
        elif record_type == 'queue':
            self._flush_queues()
        else:
            raise Exception("Unknown record type '{record_type}'".format(
                record_type=record_type))

    def _flush_all(self): self.utils.mc_db.flush_mc_db()

    def _flush_flows(self): self.utils.mc_db.delete_items(item_type='flow')

    def _flush_jobs(self): self.utils.mc_db.delete_items(item_type='job')

    def _flush_queues(self): self.utils.mc_db.delete_items(item_type='queue')


Subcommand = FlushMcRecordsSubcommand
