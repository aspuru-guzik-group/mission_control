import sqlalchemy as _sqla
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mc.utils import update_helper
from .query_builder import QueryBuilder


class Db(object):
    ITEM_TYPES = ['job', 'flow', 'queue']

    class ItemNotFoundError(Exception):
        pass

    def __init__(self, engine=None, db_uri=None, schema=None,
                 ensure_tables=False):
        if engine:
            self.engine = engine
        else:
            self.db_uri = db_uri
        self.schema = schema or self._get_default_schema()
        if ensure_tables:
            self.ensure_tables()
        self.query_builder = QueryBuilder()

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

    def ensure_tables(self): self.create_tables()

    def create_tables(self):
        self.schema.metadata.create_all(self.engine)

    def drop_tables(self):
        self.schema.metadata.drop_all(self.engine)

    @property
    def models(self):
        return self.schema.models

    def create_item(self, item_type=None, item_kwargs=None):
        Model = self.get_model_for_item_type(item_type)
        with self.session.begin_nested():
            item = Model(**item_kwargs)
            self.session.add(item)
        return self.item_to_dict(item)

    def get_model_for_item_type(self, item_type):
        return getattr(self.models, item_type.title())

    def query_items(self, item_type=None, query=None):
        """Get items of item_type that match the given query.

        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict, optional): a dict in this shape:
                ::

                    {'filters': [filter_dict_1, ...filter_dict_n]}

                where a filter_dict has this shape: ::
                    {'field': 'prop_name', 'op': 'op_name',
                     'arg': 'op argument'}

        Returns:
            items (list): a list of retrieved items.
        """
        q = self.generate_item_query(item_type=item_type, query_spec=query)
        return self.items_to_dicts(items=q)

    def generate_item_query(self, item_type=None, query_spec=None):
        Model = self.get_model_for_item_type(item_type)
        base_query = self.session.query(Model)
        item_query = self.query_builder.alter_query_per_query_spec(
            query=base_query, query_spec=query_spec)
        return item_query

    def patch_item(self, item_type=None, key=None, patches=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            key (str): the item's key.
            patches (dict): a dict of item props to update.

        Returns:
            patched_item <dict>: the patched item.
        """
        Model = self.get_model_for_item_type(item_type)
        with self.session.begin_nested():
            item = self.session.query(Model).filter_by(key=key).first()
            for field, value in patches.items():
                setattr(item, field, value)
            self.session.add(item)
        return item.to_dict()

    def flush(self, item_types=None):
        """clear the db tables."""
        for item_type in (item_types or self.ITEM_TYPES):
            Model = self.get_model_for_item_type(item_type)
            self.session.query(Model).delete()

    def items_to_dicts(self, items):
        return [self.item_to_dict(item) for item in items]

    def item_to_dict(self, item): return item.to_dict()

    def get_lock_count_subquery(self):
        lock_query = (
            self.session.query(
                _sqla.func.count(self.models.Lock.key).label('lock_count'),
                self.models.Lock.lockee_key.label('lockee_key')
            )
            .group_by(self.models.Lock.lockee_key)
        )
        return lock_query.subquery()

    def delete_items(self, item_type=None, query=None):
        """
        Args:
            item_type (str): one of :attr:`.ITEM_TYPES`
            query (dict): a query dict

        Returns:
            None
        """
        q = self.generate_item_query(item_type=item_type, query_spec=query)
        return {'num_deleted': q.delete(synchronize_session='fetch')}

    def get_item_by_key(self, item_type=None, key=None):
        """
        Args:
            item_type (str): item_type <str>: one of :attr:`.ITEM_TYPES`
            key (str): the item's key

        Raises:
            ItemNotFoundError
        """
        try:
            return self.query_items(item_type=item_type, query={
                'filters': [{'field': 'key', 'op': '=', 'arg': key}]
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
        if items_to_claim:
            claimed_items = self.patch_items(
                item_type=queue['queue_spec']['item_type'],
                keyed_patches={
                    item['key']: {'claimed': True} for item in items_to_claim
                }
            ).values()
        else:
            claimed_items = []
        return {'items': claimed_items}

    def get_queue_items_to_claim(self, queue=None):
        """
        Args:
            queue (dic): a queue record

        Returns:
            items (list): a list of items that match the queue's query.
        """
        queue_item_type = queue['queue_spec']['item_type']
        if queue_item_type == 'flow':
            claim_fn = self.get_flow_queue_items_to_claim
        else:
            claim_fn = self.default_get_queue_items_to_claim
        return claim_fn(queue=queue)

    def get_flow_queue_items_to_claim(self, queue=None):
        """Gets flow queue items.

        Checks for lock records on items.
        """
        Flow = self.models.Flow
        query = self.session.query(Flow)
        query = self.query_builder.alter_query_per_query_spec(
            query=query,
            query_spec={'filters': self.get_default_claiming_filters()}
        )
        lock_count_subquery = self.get_lock_count_subquery()
        query = (
            query.join(
                lock_count_subquery,
                (Flow.key == lock_count_subquery.c.lockee_key),
                isouter=True,
            )
            .filter(
                (Flow.num_tickable_tasks.is_(None))
                | (lock_count_subquery.c.lock_count.is_(None))
                | (Flow.num_tickable_tasks > lock_count_subquery.c.lock_count)
            )
        )
        return self.items_to_dicts(items=query)

    def get_default_claiming_filters(self):
        """
        Returns:
            filters (list): a list of default filters to use for claiming queue
                items.
        """
        return [
            {'field': 'claimed', 'op': '=', 'arg': False},
            {'field': 'status', 'op': 'IN', 'arg': ['PENDING', 'RUNNING']}
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
        return self.query_items(
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
            item_kwargs={'lockee_key': lockee_key, 'locker_key': locker_key}
        )

    def release_locks(self, locker_keys=None):
        """Release locks.

        Args:
            locker_key (str): key for the item that holds the lock.

        Returns:
            None
        """
        return self.delete_items(
            item_type='lock',
            query={
                'filters': [
                    {'field': 'locker_key', 'op': 'IN', 'arg': locker_keys},
                ]
            }
        )

    def upsert(self, key=None, updates=None, model_type=None, commit=True):
        model_type = key.split(':')[0].title()
        model = getattr(self.models, model_type)
        instance, created = self.get_or_create_instance(key=key, model=model)
        updates = self._alter_updates(updates)
        update_helper.update(instance, updates)
        self.session.merge(instance)
        if commit:
            self.session.commit()

    def _alter_updates(self, updates):
        return [self._alter_update(update) for update in updates]

    def _alter_update(self, update):
        calls_that_need_session = ['add_parents_by_key',
                                   'add_ancestors_by_key']
        if (
            update[1] == '$call' and
            update[0].split('.')[-1] in calls_that_need_session
        ):
            altered_update = self._alter_update_that_needs_session(update)
        else:
            altered_update = update
        return altered_update

    def _alter_update_that_needs_session(self, update):
        session_kwargs = {'session': self.session}
        params = [*update[2:]]
        if len(params) == 0:
            params = [[], session_kwargs]
        elif len(params) == 1:
            params.append(session_kwargs)
        elif len(params) == 2:
            params[1] = {**params[1], **session_kwargs}
        altered_update = [update[0], update[1], *params]
        return altered_update

    def get_or_create_instance(self, key=None, model=None):
        instance = self.session.query(model).filter_by(key=key).first()
        if instance:
            return instance, False
        else:
            model_kwargs = {'key': key}
            if model is self.models.Ent:
                _, ent_type, key_body = key.split(':', maxsplit=2)
                model_kwargs = {'key': key, 'ent_type': ent_type}
            instance = model(**model_kwargs)
            self.session.add(instance)
            return instance, True

    def execute_actions(self, actions=None):
        results = []
        with self.session.begin_nested():
            for action in actions or []:
                result = self.execute_action(action=action, commit=False)
                results.append(result)
        self.session.commit()
        return results

    def execute_action(self, action=None, commit=False):
        params = action.get('params', {})
        if action['type'] == 'upsert':
            fn = self.upsert
            params = {**params, 'commit': commit}
        return fn(**params)
