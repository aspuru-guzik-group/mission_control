import json

class BaseDao(object):
    ITEM_TYPES = ['Job', 'Flow', 'Queue']

    class ItemNotFoundError(Exception): pass

    def get_items(self, item_type=None, query=None): raise NotImplementedError

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

    def get_item_by_key(self, item_type=None, key=None):
        try:
            return self.get_items(item_type=item_type, query={
                'filters': [
                    {'field': 'key', 'operator': '=', 'value': key}
                ]
            })[0]
        except IndexError as exc:
            error_details = "item_type '{item_type}', key '{key}'".format(
                item_type=item_type, key=key)
            raise self.ItemNotFoundError(error_details)

    def patch_items(self, item_type=None, keyed_patches=None):
        return {
            key: self.patch_item(item_type=item_type, key=key, patches=patches)
            for key, patches in keyed_patches.items()
        }

    def deserialize_value(self, serialized_value=None):
        if serialized_value is None: return None
        return json.loads(serialized_value)

    def serialize_value(self, value=None): return json.dumps(value)

    def claim_queue_items(self, queue_key=None, **kwargs):
        queue = self.get_item_by_key(item_type='Queue', key=queue_key)
        items_to_claim = self.get_queue_items_to_claim(queue=queue)
        claimed_items = self.patch_items(
            item_type=queue['queue_spec']['item_type'],
            keyed_patches={item['key']: {'claimed': True}
                           for item in items_to_claim})
        return {'items': claimed_items.values()}

    def get_queue_items_to_claim(self, queue=None):
        queue_item_type = queue['queue_spec']['item_type']
        if queue_item_type == 'Flow':
            return self.get_flow_queue_items_to_claim(queue=queue)
        return self.default_get_queue_items_to_claim(queue=queue)

    def get_default_claiming_filters(self):
        return [
            {'field': 'claimed', 'operator': '=', 'value': False},
            {'field': 'status', 'operator': 'IN',
             'value': ['PENDING', 'RUNNING']}
        ]

    def get_flow_queue_items_to_claim(self, queue=None):
        raise NotImplementedError

    def default_get_queue_items_to_claim(self, queue=None, filters=None):
        return self.get_items(
            item_type=queue['queue_spec']['item_type'],
            query={'filters': self.get_default_claiming_filters()}
        )
