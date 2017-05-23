import logging
from sqlalchemy import create_engine

from .base_dao import BaseDao

class SaDao(BaseDao):
    def __init__(self, schema=None, db_uri=None, engine=None, logger=None):
        self.logger = logger or logging
        self.db_uri = db_uri
        self.schema = schema or self.get_default_schema()
        if engine: self.engine = engine

    def get_default_schema(self):
        from mc.orm import sqlalchemy as _mc_sa
        return _mc_sa.generate_schema()

    def create_tables(self):
        self.schema['metadata'].create_all(self.engine)

    def create_item(self, item_type=None, kwargs=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.insert().values(kwargs).return_defaults()
        key = self.engine.execute(statement).inserted_primary_key[0]
        return {
            **kwargs,
            **self.get_item_by_key(item_type=item_type, key=key)
        }

    def get_item_table(self, item_type=None):
        return self.schema['tables'][item_type]

    @property
    def engine(self):
        if not hasattr(self, '_engine'):
            self._engine = create_engine(self.db_uri)
        return self._engine

    @engine.setter
    def engine(self, value):
        self._engine = value

    def get_items(self, item_type=None, query=None):
        statement = self.get_items_statement(item_type=item_type, query=query)
        return [self.row_to_dict(row)
                for row in self.engine.execute(statement).fetchall()]
    
    def row_to_dict(self, row): return dict(zip(row.keys(), row))

    def get_items_statement(self, item_type=None, query=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.select()
        statement = self.add_wheres_to_statement(query=query,
                                                 statement=statement,
                                                 columns=table.columns)
        return statement

    def add_wheres_to_statement(self, query=None, statement=None,
                                columns=None):
        for _filter in query.get('filters', []):
            col = columns[_filter['field']]
            clause = col.op(_filter['operator'])(_filter['value'])
            statement = statement.where(clause)
        return statement

    def patch_item(self, item_type=None, key=None, patches=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.update().values(patches)\
                .where(table.columns.key == key)
        self.engine.execute(statement)
        return self.get_item_by_key(item_type=item_type, key=key)

    def flush_mc_db(self, item_types=None):
        for item_type in (item_types or self.ITEM_TYPES):
            table = self.get_item_table(item_type=item_type)
            self.engine.execute(table.delete())

    def get_queue_spec_for_queue_key(self, queue_key=None):
        queue = self.get_item_by_key(item_type='Queue', key=queue_key)
        return self.deserialize_value(queue['queue_spec'])

    def get_queue_items(self, queue_key=None, query=None):
        self.fail()
