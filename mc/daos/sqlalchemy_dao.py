from types import SimpleNamespace
import logging
import sqlalchemy as _sqla

from .base_dao import BaseDao

class SqlAlchemyDao(BaseDao):
    """A SqlAlchemy implementation of a MissionControl DAO."""
    def __init__(self, sa_schema=None, marsh_schemas=None, db_uri=None,
                 engine=None, logger=None):
        """
        Args:
            sa_schema (dict, optional): a dict representing a sqlalchemy schema,
                in this shape: ::

                    {
                        'metadata': sqlalchemy metadata obj,
                        'tables': a dict of sqlalchemy tables objects
                    }

                If empty, will use default schema.
            marsh_schemas (dict, optional): a dict representing Marshmallow
                serializers. If empty, will use default schema.
            db_uri (str, optional): SqlAlchemy db_uri string. 'sqlite://'.
                If engine is provided, this will be ignored.
            engine (SqlAlchemy engine, optional): SqlAlchemy engine object.
                If empty, will attempt to create engine from db_uri.
            logger (logging.Logger, optional): Logger instance. If empty,
                will create logger.
        """
        self.logger = logger or logging
        self.db_uri = db_uri
        self.sa_schema = sa_schema or self.get_default_sa_schema()
        self.marsh_schemas = marsh_schemas or self.get_default_marsh_schemas()
        if engine: self.engine = engine

    def get_default_sa_schema(self):
        from mc.orm import sqlalchemy as _mc_sa
        return _mc_sa.generate_schema()

    def get_default_marsh_schemas(self):
        from mc.orm import marshmallow as _mc_marsh
        return _mc_marsh.sa_schema_to_marsh_schemas(sa_schema=self.sa_schema)

    def ensure_tables(self): self.create_tables()

    def create_tables(self):
        self.sa_schema['metadata'].create_all(self.engine)

    def drop_tables(self):
        self.sa_schema['metadata'].drop_all(self.engine)

    def create_item(self, item_type=None, item_kwargs=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.insert().values(item_kwargs).return_defaults()
        key = self.engine.execute(statement).inserted_primary_key[0]
        return self.get_item_by_key(item_type=item_type, key=key)

    def serialize_items(self, item_type=None, items=None):
        serializer = self.get_marsh_schema(item_type=item_type)
        return [serializer.dump(SimpleNamespace(**item)).data for item in items]

    def get_marsh_schema(self, item_type=None):
        return self.marsh_schemas[item_type]

    def get_item_table(self, item_type=None):
        return self.sa_schema['tables'][item_type]

    @property
    def engine(self):
        if not hasattr(self, '_engine'):
            self._engine = _sqla.create_engine(self.db_uri)
        return self._engine

    @engine.setter
    def engine(self, value):
        self._engine = value

    def get_items(self, item_type=None, query=None):
        statement = self.get_items_statement(item_type=item_type, query=query)
        return self.statement_to_dicts(statement=statement)

    def statement_to_dicts(self, statement=None):
        """execute statement and return results as dicts."""
        return [self.row_to_dict(row)
                for row in self.engine.execute(statement).fetchall()]

    def row_to_dict(self, row): return dict(zip(row.keys(), row))

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
        for _filter in query.get('filters', []):
            clause = None
            col = columns[_filter['prop']]
            op = _filter['op']
            arg = _filter['arg']
            negate = False
            if op.startswith('!'):
                op  = op[2:]
                negate = True
            if op == 'IN':
                op_fn = col.in_
                if not arg: clause = _sqla.false()
            else: op_fn = col.op(op)
            if clause is None: clause = op_fn(arg)
            if negate: clause = ~clause
            statement = statement.where(clause)
        return statement

    def patch_item(self, item_type=None, key=None, patches=None):
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
                (joined.c.num_tickable_tasks == None)
                | (joined.c.lock_count == None)
                | (joined.c.num_tickable_tasks > joined.c.lock_count)
            )
        )
        return self.statement_to_dicts(statement=statement)

    def get_lock_count_subq(self):
        lock_table = self.get_item_table(item_type='Lock')
        subq = (
            _sqla.select([_sqla.func.count('key').label('lock_count'),
                          lock_table.c.lockee_key])
            .group_by(lock_table.c.lockee_key)
        )
        return subq

    def delete_items(self, item_type=None, query=None):
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
