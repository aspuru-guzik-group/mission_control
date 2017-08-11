def dump_flows(self, args=None, kwargs=None, unparsed_args=None):
    if 'keys_to_exclude' not in kwargs:
        kwargs = {**kwargs, 'keys_to_exclude': {'graph'}}
    self._dump_items(item_type='flow', kwargs=kwargs)


def dump_jobs(self, args=None, kwargs=None, unparsed_args=None):
    if 'keys_to_exclude' not in kwargs:
        kwargs = {**kwargs, 'keys_to_exclude': {'data'}}
    self._dump_items(item_type='job', **kwargs)


def _dump_items(self, item_type=None, keys_to_exclude=None, filters=None,
                **kwargs):
    print('==== ' + item_type.upper() + ' ====')
    mc_db = self._get_mc_db()
    keys_to_exclude = keys_to_exclude or {}
    for item in mc_db.get_items(item_type=item_type):
        if not all([filter_(item) for filter_ in (filters or [])]):
            continue
        self._dump_item(item=item, keys_to_exclude=keys_to_exclude)
        print('-' * 10)


def _dump_item(self, item=None, keys_to_exclude=None):
    for key, value in item.items():
        if key not in keys_to_exclude:
            print("{key}: {value}".format(key=key, value=value))


def dump_locks(self, args=None, kwargs=None, unparsed_args=None):
    self._dump_items(item_type='lock', **kwargs)
