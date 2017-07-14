def ensure_queues(self, args=None, kwargs=None, unparsed_args=None):
    mc_dao = self._get_mc_dao()
    for item_type in ['flow', 'job']:
        queue_key = self.cfg.get(item_type.upper() + '_QUEUE_KEY')
        try:
            mc_dao.create_item(
                item_type='queue',
                item_kwargs={
                    'key': queue_key,
                    'queue_spec': {'item_type': item_type},
                }
            )
            self.logger.info("Created {item_type} queue".format(
                item_type=item_type))
        except mc_dao.IntegrityError:
            self.logger.info(("Queue with key '{queue_key}' already"
                              " exists.").format(queue_key=queue_key))
