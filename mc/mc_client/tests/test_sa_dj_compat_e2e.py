import tempfile

from django.db import connections
import django.test

from .. import sa_dao
from .. import dj_dao


class BaseTestCase(django.test.TestCase):
    def setUp(self):
        self.db_cfg = {
            'id': 'tmp',
            'NAME': tempfile.mkstemp()[1],
            'ENGINE': 'django.db.backends.sqlite3',
        }
        connections.databases[self.db_cfg['id']] = self.db_cfg
        self.dj_dao = dj_dao.DjDao(db_id=self.db_cfg['id'])
        db_uri = 'sqlite:///{}'.format(self.db_cfg['NAME'])
        self.sa_dao = sa_dao.SaDao(db_uri=db_uri)
        self.sa_dao.create_tables()

    def tearDown(self):
        del connections.databases[self.db_cfg['id']]

class FlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_kwargs = {'label': 'initial_label'}

    def _create_flow(self, dao=None):
        return dao.create_item(item_type='Flow', kwargs=self.flow_kwargs)

    def test_dj_create_sa_read(self):
        flow_from_dj = self._create_flow(dao=self.dj_dao)
        flow_from_sa = self.sa_dao.get_item_by_key(item_type='Flow',
                                                   key=flow_from_dj['key'])
        print("dj: ", flow_from_dj)
        print("sa: ", flow_from_sa)
        self.assertEqual(flow_from_dj, flow_from_sa)
