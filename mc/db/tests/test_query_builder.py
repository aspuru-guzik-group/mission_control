import uuid
import unittest

import sqlalchemy as _sqla
import sqlalchemy.orm as _orm
from sqlalchemy.ext.declarative import declarative_base

from .. import query_builder


Base = declarative_base()


class MyModel(Base):
    __tablename__ = 'my_model'
    key = _sqla.Column(_sqla.String, primary_key=True,
                       default=(lambda ctx: str(uuid.uuid4())))
    int_attr = _sqla.Column(_sqla.Integer)
    str_attr = _sqla.Column(_sqla.String)


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = _sqla.create_engine('sqlite://')
        Base.metadata.create_all(self.engine)
        self.session = _orm.Session(self.engine)
        self.query_builder = query_builder.QueryBuilder()

    def generate_instance(self, **kwargs):
        return MyModel(**kwargs)


class FiltersTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.instances = {}
        for i in range(3):
            instance = self.generate_instance(
                key='instance_%s' % i,
                int_attr=i,
                str_attr=str(i)
            )
            self.instances[instance.key] = instance
        self.session.add_all(list(self.instances.values()))
        self.session.commit()

    def test_filters(self):
        q = self.session.query(MyModel)
        q = self.query_builder.alter_query_per_filters(
            query=q,
            filters=[
                {'field': 'int_attr', 'op': '=', 'arg': 1},
                {'field': 'str_attr', 'op': '=', 'arg': '1'},
                {'field': 'str_attr', 'op': '!=', 'arg': 'foo'},
            ]
        )
        results = q.all()
        self.assertEqual(set(results), set([self.instances['instance_1']]))
