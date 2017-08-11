from ._base_houston_subcommand import BaseHoustonSubcommand


class EnsureQueuesSubcommand(BaseHoustonSubcommand):
    def _run(self):
        self.utils.ensure_queues()
        for item_type in ['flow', 'job']:
            queue_cfg = self.utils.cfg.get(item_type.upper() + '_QUEUE')
            queue_key = queue_cfg['key']
            try:
                self.utils.mc_db.create_item(
                    item_type='queue',
                    item_kwargs={
                        'key': queue_key,
                        **queue_cfg.get('queue_kwargs', {})
                    }
                )
                self.logger.info("Created {item_type} queue".format(
                    item_type=item_type))
            except self.utils.mc_db.IntegrityError:
                self.logger.info(("Queue with key '{queue_key}' already"
                                  " exists.").format(queue_key=queue_key))


Subcommand = EnsureQueuesSubcommand
