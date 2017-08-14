import unittest

import sqlalchemy as _sqla
import sqlalchemy.orm as _sqla_orm

from .. import models


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = _sqla.create_engine('sqlite://')
        models.utils.Base.metadata.create_all(self.engine)
        self.Session = _sqla_orm.sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.request_type = 'some_request_type'
        self.request_tag = 'some_request_tag'
        self.instance_key = 'some_instance_key'
        self.params = {'param_%s' % i: 'value_%s' % i for i in range(3)}


class CreateRequestTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.request = models.Request(
            request_type=self.request_type,
            request_tag=self.request_tag,
            instance_key=self.instance_key,
            params=self.params
        )
        self.session.add(self.request)
        self.session.commit()

    def test_has_timestamps(self):
        self.assertTrue(self.request.created is not None)
        self.assertTrue(self.request.modified is not None)

    def test_key_startswith_request(self):
        self.assertTrue(self.request.key.startswith('request:'))

    def test_has_request_type(self):
        self.assertTrue(self.request.request_type, self.request_type)

    def test_has_request_tag(self):
        self.assertTrue(self.request.request_tag, self.request_tag)

    def test_has_instance_key(self):
        self.assertTrue(self.request.instance_key, self.instance_key)

    def test_has_params(self):
        self.assertTrue(self.request.params, self.params)

    def test_has_default_status(self):
        self.assertTrue(self.request.status, 'PENDING')
