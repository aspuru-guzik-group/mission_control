from types import SimpleNamespace
import logging
from sqlalchemy import create_engine

from .base_dao import BaseDao

class SaDao(BaseDao):
    def __init__(self, sa_schema=None, marsh_schemas=None, db_uri=None,
                 engine=None, logger=None):
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

    def create_tables(self):
        self.sa_schema['metadata'].create_all(self.engine)

    def drop_tables(self):
        self.sa_schema['metadata'].drop_all(self.engine)

    def create_item(self, item_type=None, kwargs=None):
        table = self.get_item_table(item_type=item_type)
        statement = table.insert().values(kwargs).return_defaults()
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
            self._engine = create_engine(self.db_uri)
        return self._engine

    @engine.setter
    def engine(self, value):
        self._engine = value

    def get_items(self, item_type=None, query=None):
        statement = self.get_items_statement(item_type=item_type, query=query)
        items = [self.row_to_dict(row)
                 for row in self.engine.execute(statement).fetchall()]
        return self.serialize_items(item_type=item_type, items=items)
    
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
            operator = _filter['operator']
            if operator == 'IN': op_fn = col.in_
            else: op_fn = col.op(operator)
            clause = op_fn(_filter['value'])
            statement = statement.where(clause)
        return statement

    def patch_item(self, item_type=None, key=None, patches=None):
        table = self.get_item_table(item_type=item_type)
        statement = \
                table.update().where(table.columns.key == key).values(patches)
        self.engine.execute(statement)
        return self.get_item_by_key(item_type=item_type, key=key)

    def flush_mc_db(self, item_types=None):
        for item_type in (item_types or self.ITEM_TYPES):
            table = self.get_item_table(item_type=item_type)
            self.engine.execute(table.delete())

    def get_queue_spec_for_queue_key(self, queue_key=None):
        queue = self.get_item_by_key(item_type='Queue', key=queue_key)
        return self.deserialize_value(queue['queue_spec'])
