import sqlalchemy as _sqla
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


class Db(object):
    ITEM_TYPES = ['job', 'flow', 'queue']

    class InvalidQueryError(Exception):
        pass

    class InvalidFilterError(Exception):
        pass

    class ItemNotFoundError(Exception):
        pass

    class IntegrityError(Exception):
        pass

    def __init__(self, engine=None, db_uri=None, schema=None, ensure=False):
        if engine:
            self.engine = engine
        else:
            self.db_uri = db_uri
        self.schema = schema or self._get_default_schema()
        if ensure:
            self.ensure_db()

    def _get_default_schema(self):
        from . import schema
        return schema

    @property
    def engine(self):
        if not hasattr(self, '_engine'):
            db_uri = self.db_uri
            if callable(db_uri):
                db_uri = db_uri()
            self._engine = create_engine(db_uri)
        return self._engine

    @engine.setter
    def engine(self, value): self._engine = value

    @property
    def session(self):
        if not hasattr(self, '_session'):
            self._session = self.Session()
        return self._session

    @session.setter
    def session(self, value): self._session = value

    @property
    def Session(self):
        if not hasattr(self, '_Session'):
            self._Session = sessionmaker(bind=self.engine)
        return self._Session

    @Session.setter
    def Session(self, value): self._Session = value

    def ensure_db(self):
        self.schema.metadata.create_all(self.engine)

    @property
    def models(self):
        return self.schema.models

    def get_or_create_instance(self, key=None, model=None):
        instance = self.session.query(model).filter_by(key=key).first()
        if instance:
            return instance, False
        else:
            instance = model(key=key)
            self.session.add(instance)
            return instance, True

    def create_item(self, item_type=None, item_kwargs=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.insert().values(item_kwargs).return_defaults()
        try:
            key = self.engine.execute(statement).inserted_primary_key[0]
        except _sqla.exc.IntegrityError as exc:
            raise self.IntegrityError() from exc
        return self.get_item_by_key(item_type=item_type, key=key)

    def get_items(self, item_type=None, query=None):
        """Get items of item_type that match the given query.

        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict, optional): a dict in this shape:
                ::

                    {'filters': [filter_dict_1, ...filter_dict_n]}

                where a filter_dict has this shape: ::
                    {'prop': 'prop_name', 'op': 'op_name',
                     'arg': 'op argument'}

        Returns:
            items (list): a list of retrieved items.
        """
        q = self.get_items_statement(item_type=item_type, query=query)
        return self.statement_to_dicts(statement=statement)

    def get_items_statement(self, item_type=None, query=None):
        """get statement for items query."""
        query = query or {}
        table = self.get_item_table(item_type=item_type)
        statement = table.select()
        statement = self.add_wheres_to_statement(query=query,
                                                 statement=statement,
                                                 columns=table.columns)
        return statement

    def add_wheres_to_statement(self, query=None, statement=None,
                                columns=None):
        for filter_ in query.get('filters', []):
            clause = None
            col = columns[filter_['prop']]
            op = filter_['op']
            arg = filter_['arg']
            negate = False
            if op.startswith('!'):
                op = op[2:]
                negate = True
            if op == 'IN':
                op_fn = col.in_
                if not arg:
                    clause = _sqla.false()
            else:
                op_fn = col.op(op)
            if clause is None:
                clause = op_fn(arg)
            if negate:
                clause = ~clause
            statement = statement.where(clause)
        return statement

    def patch_item(self, item_type=None, key=None, patches=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            key (str): the item's key.
            patches (dict): a dict of item props to update.

        Returns:
            patched_item <dict>: the patched item.
        """
        table = self.get_item_table(item_type=item_type)
        statement = (
            table.update()
            .where(table.columns.key == key)
            .values(patches)
        )
        self.engine.execute(statement)
        return self.get_item_by_key(item_type=item_type, key=key)

    def flush_mc_db(self, item_types=None):
        """clear the db tables."""
        for item_type in (item_types or self.ITEM_TYPES):
            table = self.get_item_table(item_type=item_type)
            self.engine.execute(table.delete())

    def get_flow_queue_items_to_claim(self, queue=None):
        """Gets flow queue items.

        Checks for lock records on items.
        """
        filtered_flows_q = self.get_items_statement(
            item_type=queue['queue_spec']['item_type'],
            query={'filters': self.get_default_claiming_filters()}
        )
        lock_count_subq = self.get_lock_count_subq()
        joined = filtered_flows_q.join(
            lock_count_subq,
            isouter=True,
            onclause=(filtered_flows_q.c.key == lock_count_subq.c.lockee_key)
        )
        statement = (
            _sqla.select('*')
            .select_from(joined)
            .where(
                (joined.c.num_tickable_tasks == None)  # noqa
                | (joined.c.lock_count == None)
                | (joined.c.num_tickable_tasks > joined.c.lock_count)
            )
        )
        return self.statement_to_dicts(statement=statement)

    def get_lock_count_subq(self):
        lock_table = self.get_item_table(item_type='lock')
        subq = (
            _sqla.select([_sqla.func.count('key').label('lock_count'),
                          lock_table.c.lockee_key])
            .group_by(lock_table.c.lockee_key)
        )
        return subq

    def delete_items(self, item_type=None, query=None):
        """
        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict): a query dict

        Returns:
            None
        """
        item_table = self.get_item_table(item_type=item_type)
        keys_subq = (
            self.get_items_statement(item_type=item_type, query=query)
            .with_only_columns([item_table.c.key])
        )
        delete_statement = (
            item_table.delete()
            .where(item_table.c.key.in_(keys_subq))
        )
        self.engine.execute(delete_statement)

    def validate_query(self, query=None):
        """
        Args:
            query (dict): a query as per get_items.

        Raises:
            InvalidQueryError
        """
        filter_errors = []
        for filter_ in query.get('filters' or []):
            try:
                self.validate_query_filter(filter_=filter_)
            except self.InvalidFilterError as exc:
                filter_errors.append(exc)
        if filter_errors:
            error = "Invalid filters: {filter_errors}".format(
                filter_errors="\n".join([str(filter_error)
                                         for filter_error in filter_errors])
            )
            raise self.InvalidQueryError(error)

    def validate_query_filter(self, filter_=None):
        valid_ops = ['=', 'IN']
        if filter_['op'] not in valid_ops:
            error = ("Invalid op '{op}',  valid ops are: {valid_ops}").format(
                op=filter_['op'], valid_ops=valid_ops)
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
            patched_items (dict): a dictionary of patched results, keyed by
                item keys
        """
        return {
            key: self.patch_item(item_type=item_type, key=key, patches=patches)
            for key, patches in keyed_patches.items()
        }

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
        queue = self.get_item_by_key(item_type='queue', key=queue_key)
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
        if queue_item_type == 'flow':
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
            item_type='lock',
            item_kwargs={'lockee_key': lockee_key, 'locker_key': locker_key}
        )

    def release_locks(self, locker_keys=None):
        """Release locks.

        Args:
            locker_key (str): key for the item that holds the lock.

        Returns:
            None
        """
        self.delete_items(
            item_type='lock',
            query={
                'filters': [
                    {'prop': 'locker_key', 'op': 'IN', 'arg': locker_keys},
                ]
            }
        )
