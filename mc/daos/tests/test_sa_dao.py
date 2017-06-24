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

    def setup_module_mocks(self, attrs=None, module=sa_dao):
        patchers = {attr: patch.object(module, attr) for attr in attrs}
        mocks = {}
        for attr, patcher in patchers.items():
            self.addCleanup(patcher.stop)
            mocks[attr] = patcher.start()
        return mocks

    def setup_dao_mocks(self, attrs=None, mock_factory=MagicMock):
        for attr in attrs: setattr(self.dao, attr, mock_factory())

    def generate_mocks(self, n=3): return [MagicMock() for i in range(n)]

class CreateItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_dao_mocks(attrs=['get_item_table', 'get_item_by_key',
                                    'execute'])
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

    @patch.object(sa_dao, '_sqla')
    def test_creates_engine_if_not_set(self, _sqla):
        del self.dao._engine
        self.assertEqual(self.dao.engine, _sqla.create_engine.return_value)

    def test_sets_engine(self):
        value = MagicMock()
        self.dao.engine = value
        self.assertEqual(self.dao._engine, value)

class GetItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_dao_mocks(attrs=['get_items_statement',
                                    'statement_to_dicts'])
        self.result = self.dao.get_items(item_type=self.item_type,
                                         query=self.query)

    def test_gets_items_statement(self):
        self.assertEqual(self.dao.get_items_statement.call_args,
                         call(item_type=self.item_type, query=self.query))

    def test_returns_results(self):
        self.assertEqual(
            self.dao.statement_to_dicts.call_args,
            call(statement=self.dao.get_items_statement.return_value)
        )
        self.assertEqual(self.result, self.dao.statement_to_dicts.return_value)

class GetItemsStatementTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_dao_mocks(attrs=['get_item_table',
                                    'add_wheres_to_statement'])
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

class StatementToDictsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.statement = MagicMock()
        self.setup_dao_mocks(attrs=['row_to-dict'])
        self.mock_rows = [MagicMock() for i in range(3)]
        self.dao.engine.execute.return_value.fetchall.return_value = \
                self.mock_rows
        self.result = self.dao.statement_to_dicts(statement=self.statement)

    def test_executes_statement(self):
        self.assertEqual(self.dao.engine.execute.call_args,
                         call(self.statement))

    def test_returns_rows_as_dicts(self):
        expected_items = [
            self.dao.row_to_dict(row)
            for row in self.dao.engine.execute.return_value.fetchall()
        ]
        self.assertEqual(self.result, expected_items)

class PatchItemTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = MagicMock()
        self.patches = MagicMock()
        self.setup_dao_mocks(attrs=['get_item_table', 'get_item_by_key'])
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
        self.setup_dao_mocks(attrs=['get_item_table'])
        self.dao.flush_mc_db()

    def test_calls_delete_for_each_item_table(self):
        for i, item_type in enumerate(self.dao.ITEM_TYPES):
            self.assertEqual(self.dao.get_item_table.call_args_list[i],
                             call(item_type=item_type))
            expected_table = self.dao.get_item_table.return_value
            self.assertEqual(self.dao.engine.execute.call_args,
                             call(expected_table.delete()))

class DeleteItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_dao_mocks(attrs=['get_item_table', 'get_items_statement'])
        self.item_type = MagicMock()
        self.query = MagicMock()
        self.result = self.dao.delete_items(item_type=self.item_type,
                                            query=self.query)

    def test_gets_item_table(self):
        self.assertEqual(self.dao.get_item_table.call_args,
                         call(item_type=self.item_type))

    def test_gets_items_statement(self):
        self.assertEqual(self.dao.get_items_statement.call_args,
                         call(item_type=self.item_type, query=self.query))

    def test_builds_and_executes_delete_statement(self):
        expected_item_table = self.dao.get_item_table.return_value
        expected_items_statement = self.dao.get_items_statement.return_value
        self.assertEqual(expected_items_statement.with_only_columns.call_args,
                         call([expected_item_table.c.key]))
        expected_keys_subq = (expected_items_statement.with_only_columns
                              .return_value)
        self.assertEqual(
            expected_item_table.delete.return_value.where.call_args,
            call(expected_item_table.c.key.in_(expected_keys_subq))
        )
        expected_delete_statement = (expected_item_table.delete.return_value
                                     .where.return_value)
        self.assertEqual(self.dao.engine.execute.call_args,
                         call(expected_delete_statement))

class GetFlowQueueItemsToClaim(BaseTestCase):
    def setUp(self):
        super().setUp()

        class MagicMockWithGt(MagicMock):
            def __init__(_self, *args, **kwargs):
                super(MagicMock, _self).__init__(*args, **kwargs)
                def mock_gt(*args, **kwargs): return _self.cls()
                _self.__gt__ = mock_gt

        mock_factory = MagicMockWithGt
        self.setup_dao_mocks(attrs=['get_items_statement',
                                    'get_default_claiming_filters',
                                    'get_lock_count_subq',
                                    'statement_to_dicts'],
                             mock_factory=mock_factory)
        self.mocks = self.setup_module_mocks(attrs=['_sqla'])
        self.queue = mock_factory()
        self.expected_items_statement = \
                self.dao.get_items_statement.return_value
        self.expected_lock_subq = self.dao.get_lock_count_subq.return_value
        self.result = self.dao.get_flow_queue_items_to_claim(queue=self.queue)

    def test_gets_claiming_filters(self):
        self.assertEqual(self.dao.get_default_claiming_filters.call_args,
                         call())

    def test_gets_items_statement(self):
        self.assertEqual(
            self.dao.get_items_statement.call_args,
            call(item_type=self.queue['queue_spec']['item_type'],
                 query={'filters': self.dao.get_default_claiming_filters()}
                )
        )

    def test_gets_lock_count_subq(self):
        self.assertEqual(self.dao.get_lock_count_subq.call_args, call())


    def test_joins_items_statement_and_lock_subq(self):
        self.assertEqual(
            self.expected_items_statement.join.call_args,
            call(self.expected_lock_subq,
                 isouter=True,
                 onclause=(self.expected_items_statement.c.key == \
                           self.expected_lock_subq.c.lockee_key)
                )
        )

    def test_selects_from_join_with_where_conditions(self):
        expected_join = self.expected_items_statement.join.return_value
        self.assertEqual(self.mocks['_sqla'].select.call_args, call('*'))
        expected_select = self.mocks['_sqla'].select.return_value
        self.assertEqual(expected_select.select_from.call_args,
                         call(expected_join))
        self.assertEqual(
            expected_select.select_from.return_value.where.call_args,
            call(
                (expected_join.c.num_tickable_tasks == None)
                | (expected_join.c.lock_count == None)
                | (expected_join.c.num_tickable_tasks > \
                   expected_join.c.lock_count)
            )
        )

    def test_returns_dicts(self):
        self.assertEqual(
            self.dao.statement_to_dicts.call_args,
            call(statement=(self.mocks['_sqla'].select.return_value.select_from
                            .return_value.where.return_value))
        )
        self.assertEqual(self.result, self.dao.statement_to_dicts.return_value)

if __name__ == '__main__':
    unittest.main()


