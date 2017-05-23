import unittest
from unittest.mock import call, MagicMock, patch

from .. import sa_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.item_type = MagicMock()
        self.query = MagicMock()
        self.dao = sa_dao.SaDao(
            db_uri=MagicMock(),
            sa_schema=MagicMock(),
            marsh_schemas=MagicMock(),
            engine=MagicMock(),
        )

    def generate_mocks(self, n=3): return [MagicMock() for i in range(n)]

class CreateItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_table = MagicMock()
        self.dao.get_item_by_key = MagicMock()
        self.dao.execute = MagicMock()
        self.expected_table = self.dao.get_item_table.return_value
        self.kwargs = MagicMock()
        self.result = self.dao.create_item(item_type=self.item_type,
                                           kwargs=self.kwargs)

    def test_gets_item_table(self):
        self.assertEqual(self.dao.get_item_table.call_args,
                         call(item_type=self.item_type))

    def test_creates_insert_statement(self):
        self.assertEqual(
            self.expected_table.insert.return_value.values.call_args,
            call(self.kwargs)
        )

    def test_executes_insert_statement(self):
        expected_statement = self.expected_table.insert.return_value\
                .values.return_value.return_defaults.return_value
        self.assertEqual(self.dao.engine.execute.call_args,
                         call(expected_statement))

    def test_returns_fetched_item(self):
        expected_key = \
                self.dao.engine.execute.return_value.inserted_primary_key[0]
        self.assertEqual(self.dao.get_item_by_key.call_args,
                         call(item_type=self.item_type, key=expected_key))
        self.assertEqual(self.result, self.dao.get_item_by_key.return_value)

class GetItemTableTestCase(BaseTestCase):
    def test_gets_table_from_dict(self):
        result = self.dao.get_item_table(item_type=self.item_type)
        self.assertEqual(result, self.dao.sa_schema['tables'][self.item_type])

class EngineTestCase(BaseTestCase):
    def test_gets_engine_if_set(self):
        self.dao._engine = MagicMock()
        self.assertEqual(self.dao.engine, self.dao._engine)

    @patch.object(sa_dao, 'create_engine')
    def test_creates_engine_if_not_set(self, _create_engine):
        del self.dao._engine
        self.assertEqual(self.dao.engine, _create_engine.return_value)

    def test_sets_engine(self):
        value = MagicMock()
        self.dao.engine = value
        self.assertEqual(self.dao._engine, value)

class GetItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_items_statement = MagicMock()
        self.dao.serialize_items = MagicMock()
        self.result = self.dao.get_items(item_type=self.item_type,
                                         query=self.query)

    def test_gets_items_statement(self):
        self.assertEqual(self.dao.get_items_statement.call_args,
                         call(item_type=self.item_type, query=self.query))

    def test_returns_serialized_results(self):
        self.assertEqual(
            self.dao.engine.execute.call_args,
            call(self.dao.get_items_statement.return_value)
        )
        expected_items = [
            self.dao.row_to_dict(row)
            for row in self.dao.engine.execute.return_value.fetchall()
        ]
        self.assertEqual(self.dao.serialize_items.call_args,
                         call(item_type=self.item_type, items=expected_items))
        self.assertEqual(self.result, self.dao.serialize_items.return_value)

class GetItemsStatementTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_table = MagicMock()
        self.dao.add_wheres_to_statement = MagicMock()
        self.result = self.dao.get_items_statement(item_type=self.item_type,
                                                   query=self.query)

    def test_gets_item_table(self):
        self.assertEqual(self.dao.get_item_table.call_args,
                         call(item_type=self.item_type))

    def test_returns_statement_with_wheres(self):
        expected_table = self.dao.get_item_table.return_value
        self.assertEqual(
            self.dao.add_wheres_to_statement.call_args,
            call(statement=expected_table.select(),
                 query=self.query,
                 columns=expected_table.columns)
        )
        self.assertEqual(self.result,
                         self.dao.add_wheres_to_statement.return_value)

class PatchItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = MagicMock()
        self.patches = MagicMock()
        self.dao.get_item_table = MagicMock()
        self.dao.get_item_by_key = MagicMock()
        self.result = self.dao.patch_item(item_type=self.item_type,
                                          key=self.key, patches=self.patches)

    def test_gets_table(self):
        self.assertEqual(self.dao.get_item_table.call_args,
                         call(item_type=self.item_type))

    def test_updates_item_row(self):
        expected_table = self.dao.get_item_table.return_value
        self.assertEqual(expected_table.update().where.call_args,
                         call(expected_table.columns.key == self.key))
        self.assertEqual(
            expected_table.update().where.return_value.values.call_args,
            call(self.patches)
        )
        self.assertEqual(
            self.dao.engine.execute.call_args,
            call(expected_table.update().where.return_value.values.return_value)
        )

    def test_returns_updated_row(self):
        self.assertEqual(self.dao.get_item_by_key.call_args,
                         call(item_type=self.item_type, key=self.key))
        self.assertEqual(self.result, self.dao.get_item_by_key.return_value)

class FlushTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.get_item_table = MagicMock()
        self.dao.flush_mc_db()

    def test_calls_delete_for_each_item_table(self):
        for i, item_type in enumerate(self.dao.ITEM_TYPES):
            self.assertEqual(self.dao.get_item_table.call_args_list[i],
                             call(item_type=item_type))
            expected_table = self.dao.get_item_table.return_value
            self.assertEqual(self.dao.engine.execute.call_args,
                             call(expected_table.delete()))

if __name__ == '__main__':
    unittest.main()


