from ._base_subcommand import BaseSubcommand


class Subcommand(BaseSubcommand):
    def _run(self):
        self.utils.ensure_queues()
        for item_type in ['flow', 'job']:
            queue_cfg = self.utils.cfg.get(item_type.upper() + '_QUEUE')
            queue_key = queue_cfg['key']
            try:
                self.utils.db.create_item(
                    item_type='queue',
                    item_kwargs={
                        'key': queue_key,
                        **queue_cfg.get('queue_kwargs', {})
                    }
                )
                self.logger.info("Created {item_type} queue".format(
                    item_type=item_type))
            except self.utils.db.IntegrityError:
                self.logger.info(("Queue with key '{queue_key}' already"
                                  " exists.").format(queue_key=queue_key))
