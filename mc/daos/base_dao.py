import abc
import json

class BaseDao(abc.ABC):
    """A Base class that describes the expected DAO functions and provides
    common utilities.
    """

    ITEM_TYPES = ['Job', 'Flow', 'Queue']

    class InvalidQueryError(Exception): pass
    class InvalidFilterError(Exception): pass
    class ItemNotFoundError(Exception): pass

    @abc.abstractmethod
    def get_items(self, item_type=None, query=None):
        """Get items of item_type that match the given query.

        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict, optional): a dict in this shape:
                ::

                    {'filters': [filter_dict_1, ...filter_dict_n]}

                where a filter_dict has this shape: ::
                    {'prop': 'prop_name', 'op': 'op_name', 'arg': 'op argument'}

        Returns:
            items (list): a list of retrieved items.

        """
        raise NotImplementedError

    def validate_query(self, query=None):
        """
        Args:
            query (dict): a query as per get_items.

        Raises:
            InvalidQueryError
        """
        filter_errors = []
        for _filter in query.get('filters' or []):
            try: self.validate_query_filter(_filter=_filter)
            except self.InvalidFilterError as exc: filter_errors.append(exc)
        if filter_errors:
            error = "Invalid filters: {filter_errors}".format(
                filter_errors="\n".join([str(filter_error)
                                         for filter_error in filter_errors])
            )
            raise self.InvalidQueryError(error)

    def validate_query_filter(self, _filter=None):
        valid_ops = ['=', 'IN']
        if _filter['op'] not in valid_ops:
            error = ("Invalid op '{op}',  valid ops are: {valid_ops}").format(
                op=_filter['op'], valid_ops=valid_ops)
            self.InvalidFilterError(error)

    def get_item_by_key(self, item_type=None, key=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            key (str): the item's key

        Raises:
            ItemNotFoundError
        """
        try:
            return self.get_items(item_type=item_type, query={
                'filters': [
                    {'prop': 'key', 'op': '=', 'arg': key}
                ]
            })[0]
        except IndexError as exc:
            error_details = "item_type '{item_type}', key '{key}'".format(
                item_type=item_type, key=key)
            raise self.ItemNotFoundError(error_details)

    def patch_items(self, item_type=None, keyed_patches=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            keyed_patches (dict): a dictionary in which the keys are item_keys
                and the values are dicts of item props to update.

        Returns:
            patched_items (dict): a dictionary of patched results, keyed by item
                keys
        """
        return {
            key: self.patch_item(item_type=item_type, key=key, patches=patches)
            for key, patches in keyed_patches.items()
        }

    @abc.abstractmethod
    def patch_item(self, item_type=None, key=None, patches=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            key (str): the item's key.
            patches (dict): a dict of item props to update.

        Returns:
            patched_item <dict>: the patched item.
        """
        raise NotImplementedError

    def deserialize_value(self, serialized_value=None):
        if serialized_value is None: return None
        return json.loads(serialized_value)

    def serialize_value(self, value=None): return json.dumps(value)

    def claim_queue_items(self, queue_key=None, **kwargs):
        """
        Builds query for queue by examining queue's queue_spec.

        Args:
            queue_key (str): the queue's key

        Returns:
            claimed_items (dict): a dict of claim result, in this shape:
                ::
                    {items: [claimed_item_1, ..., claimed_item_n]}
        """
        queue = self.get_item_by_key(item_type='Queue', key=queue_key)
        items_to_claim = self.get_queue_items_to_claim(queue=queue)
        claimed_items = self.patch_items(
            item_type=queue['queue_spec']['item_type'],
            keyed_patches={item['key']: {'claimed': True}
                           for item in items_to_claim})
        return {'items': claimed_items.values()}

    def get_queue_items_to_claim(self, queue=None):
        """
        Args:
            queue (dic): a queue record

        Returns:
            items (list): a list of items that match the queue's query.
        """
        queue_item_type = queue['queue_spec']['item_type']
        if queue_item_type == 'Flow':
            return self.get_flow_queue_items_to_claim(queue=queue)
        return self.default_get_queue_items_to_claim(queue=queue)

    def get_default_claiming_filters(self):
        """
        Returns:
            filters (list): a list of default filters to use for claiming queue
                items.
        """
        return [
            {'prop': 'claimed', 'op': '=', 'arg': False},
            {'prop': 'status', 'op': 'IN', 'arg': ['PENDING', 'RUNNING']}
        ]

    @abc.abstractmethod
    def get_flow_queue_items_to_claim(self, queue=None):
        """Special handler for claiming flow queue items.
        
        Often a DAO will want to exclude a flow from being claimed if it is
        locked.

        Args:
            queue (dict): a queue record

        Returns:
            items (list): a list of items that match the queue's query.
        """
        raise NotImplementedError

    def default_get_queue_items_to_claim(self, queue=None, filters=None):
        """Default handler for claiming queue items.

        Args:
            queue (dict): a queue record
            filters <list>: filters to use for getting claimable items.

        Returns:
            items (list): a list of items that match the combination of the
                filters and the queue's queue_spec.
        """
        return self.get_items(
            item_type=queue['queue_spec']['item_type'],
            query={'filters': self.get_default_claiming_filters()}
        )

    def create_lock(self, lockee_key=None, locker_key=None):
        """Create a lock record.

        Args:
            lockee_key (str): keyfor the item being locked.
            locker_key (str): key for the item that holds the lock.

        Returns:
            lock_record (dict): a lock_record
        """
        self.create_item(
            item_type='Lock',
            kwargs={'lockee_key': lockee_key, 'locker_key': locker_key}
        )

    def release_locks(self, locker_keys=None):
        """Release locks.

        Args:
            locker_key (str): key for the item that holds the lock.

        Returns:
            None
        """
        self.delete_items(
            item_type='Lock',
            query={
                'filters': [
                    {'prop': 'locker_key', 'op': 'IN', 'arg': locker_keys},
                ]
            }
        )

    @abc.abstractmethod
    def delete_items(self, item_type=None, query=None):
        """
        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict): a query dict

        Returns:
            None
        """
        raise NotImplementedError
