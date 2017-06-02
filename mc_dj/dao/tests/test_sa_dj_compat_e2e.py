import itertools
import tempfile
import uuid

from django.db import connections
import django.test

from .. import sa_dao
from .. import dj_dao


class BaseTestCase(django.test.TestCase):
    def setUp(self):
        self.db_cfg = {
            'id': str(uuid.uuid4()),
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

class ItemTestCaseMixin(object):
    item_type = None
    item_kwargs = {'label': 'initial_label'}

    def _create_item(self, dao=None):
        return dao.create_item(item_type=self.item_type,
                               kwargs=self.item_kwargs)

    def test_create_matches_read(self):
        for dao1, dao2 in itertools.permutations([self.sa_dao, self.dj_dao], 2):
            item_from_dao1 = self._create_item(dao=dao1)
            item_from_dao2 = dao2.get_item_by_key(item_type=self.item_type,
                                                  key=item_from_dao1['key'])
            self.assertEqual(item_from_dao1, item_from_dao2)

class FlowTestCase(BaseTestCase,ItemTestCaseMixin):
    item_type = 'Flow'

class JobTestCase(BaseTestCase, ItemTestCaseMixin):
    item_type = 'Job'

class QueueTestCase(BaseTestCase, ItemTestCaseMixin):
    item_type = 'Queue'
