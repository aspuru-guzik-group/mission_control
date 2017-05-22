import logging


class DjDao(object):
    ITEM_TYPES = ['Job', 'Flow', 'Queue']

    def __init__(self, db_id=None, models=None, serializers=None, logger=None):
        self.logger = logger or logging
        self.db_id = db_id or 'default'
        self.models = models
        self.serializers = serializers

    def create_item(self, item_type=None, kwargs=None):
        model_cls = self.get_item_model_cls(item_type=item_type)
        model = model_cls.objects.db_manager(self.db_id).create(**kwargs)
        return self.serialize_item_models(item_type=item_type,
                                          item_models=[model])[0]

    def get_item_model_cls(self, item_type=None):
        return self.models[item_type]

    def serialize_item_models(self, item_type=None, item_models=None):
        serializer = self.get_item_serializer(item_type=item_type)
        return [serializer(item_model).data for item_model in item_models]

    def get_item_serializer(self, item_type=None):
        return self.serializers[item_type]

    def get_items(self, item_type=None, query=None):
        item_models = self.get_item_models(item_type=item_type, query=query)
        return self.serialize_item_models(item_type=item_type,
                                          item_models=item_models)

    def get_item_models(self, item_type=None, query=None):
        item_models = []
        queryset = self.item_query_to_queryset(item_type=item_type, query=query)
        if queryset: item_models = queryset.all()
        return item_models

    def item_query_to_queryset(self, item_type=None, query=None):
        if query: self.validate_query(query=query)
        item_model_cls = self.get_item_model_cls(item_type=item_type)
        qs = item_model_cls.objects.using(self.db_id)
        return qs.filter(**self.get_dj_filter_kwargs_for_query(query=query))

    def validate_query(self, query=None):
        for _filter in query.get('filters' or []):
            self.validate_query_filter(_filter=_filter)

    def validate_query_filter(self, _filter=None):
        valid_operators = ['=', 'IN']
        if _filter['operator'] not in valid_operators:
            raise Exception(("Invalid operator '{operator}',"
                             " valid operators are: {valid_operators}").format(
                                 operator=_filter['operator'],
                                 valid_operators=valid_operators))

    def get_dj_filter_kwargs_for_query(self, query=None):
        filter_kwargs = {}
        if not query: return filter_kwargs
        for _filter in query.get('filters', []):
            filter_kwargs.update(
                self.query_filter_to_dj_filter_kwargs(_filter=_filter))
        return filter_kwargs

    def query_filter_to_dj_filter_kwargs(self, _filter=None):
        dj_filter_kwargs = {}
        kwarg_name = _filter['field']
        if _filter['operator'] != '=':
            kwarg_name += _filter['operator'].lower()
        dj_filter_kwargs[kwarg_name] = _filter['value']
        return dj_filter_kwargs

    def patch_item(self, item_type=None, key=None, patches=None):
        model = self.get_item_model_by_key(item_type=item_type, key=key)
        patched_model = self.patch_model(model=model, patches=patches)
        return self.serialize_item_models(item_type=item_type,
                                          item_models=[patched_model])[0]

    def patch_items(self, item_type=None, keyed_patches=None):
        return {
            key: self.patch_item(item_type=item_type, patches=patches)
            for key, patches in keyed_patches.items()
        }

    def get_item_model_by_key(self, item_type=None, key=None):
        return self.get_item_models(
            item_type=item_type,
            query={
                'filters': [
                    {'field': 'key', 'operator': '=', 'value': key}
                ]
            }
        )[0]

    def patch_model(self, model=None, patches=None):
        patches = patches or {}
        for k, v in patches.items(): setattr(model, k, v)
        model.save()
        return model

    def flush_mc_db(self, item_types=None):
        for item_type in (item_types or self.ITEM_TYPES):
            model_cls = self.get_item_model_cls(item_type=item_type)
            model_cls.objects.using(self.db_id).all().delete()

    def claim_queue_items(self, queue_key=None):
        queue_model = self.get_item_model_by_key(item_type='Queue',
                                                 key=queue_key)
        queue_item_type = queue_model.queue_spec['item_type']
        items_to_claim = self.get_items(
            item_type=queue_item_type,
            query={
                'filters': [
                    {'field': 'claimed', 'operator': '=', 'value': False}
                ]
            }
        )
        claimed_items = self.patch_items(
            item_type=queue_item_type,
            keyed_patches={item['uuid']: {'claimed': True}
                           for item in items_to_claim})
        return claimed_items
