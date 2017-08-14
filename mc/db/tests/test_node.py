import unittest

import sqlalchemy as _sqla
import sqlalchemy.orm as _sqla_orm

from .. import utils


class A(utils.NodeSubclassMixin, utils.Node):
    id = _sqla.Column(_sqla.Integer, primary_key=True)
    a_name = _sqla.Column(_sqla.String)
    __mapper_args__ = {
        'polymorphic_identity': 'a',
    }


class B(utils.NodeSubclassMixin, utils.Node):
    id = _sqla.Column(_sqla.Integer, primary_key=True)
    b_name = _sqla.Column(_sqla.String)
    __mapper_args__ = {
        'polymorphic_identity': 'b',
    }


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = _sqla.create_engine('sqlite://')
        self.props = {
            'str_prop': 'str',
            'int_prop': 1,
            'bool_prop': True,
            'mapping_prop': {'some': {'nested': 'prop'}},
            'sequence_prop': ['some', 'sequence']
        }
        self.tags = {'tag_%s' % i for i in range(3)}
        self.alter_schema()
        utils.Base.metadata.create_all(self.engine)

    def alter_schema(self): pass

    @property
    def session(self):
        if not hasattr(self, '_session'):
            Session = _sqla_orm.sessionmaker(bind=self.engine)
            self._session = Session()
        return self._session


class IntraClassParentChildTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A()
        self.a2 = A(parent_nodes=[self.a1])
        self.a3 = A(parent_nodes=[self.a1])
        self.a4 = A(parent_nodes=[self.a2, self.a3])
        self.session.add_all([self.a1, self.a2, self.a3, self.a4])

    def test_nodes_have_child_nodes(self):
        self.assertEqual(self.a1.child_nodes, [self.a2, self.a3])

    def test_nodes_have_parent_nodes(self):
        self.assertEqual(self.a4.parent_nodes, [self.a2, self.a3])


class InterClassParentChildTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A()
        self.b1 = B(parent_nodes=[self.a1])
        self.a2 = A(parent_nodes=[self.b1])
        self.session.add_all([self.a1, self.b1, self.a2])

    def test_nodes_have_child_nodes(self):
        self.assertEqual(self.a1.child_nodes, [self.b1])
        self.assertEqual(self.b1.child_nodes, [self.a2])

    def test_nodes_have_parent_nodes(self):
        self.assertEqual(self.a2.parent_nodes, [self.b1])


class ScionsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A()
        self.a2 = A()
        self.a3 = A(parent_nodes=[self.a1], ancestor_nodes=[self.a1])
        self.a4 = A(parent_nodes=[self.a2], ancestor_nodes=[self.a1])
        self.a5 = A(parent_nodes=[self.a3], ancestor_nodes=[self.a1])
        self.a6 = A(parent_nodes=[self.a1, self.a2],
                    ancestor_nodes=[self.a1, self.a2])
        self.session.add_all([self.a1, self.a2, self.a3, self.a4, self.a4,
                              self.a5, self.a6])

    def test_ancestor_nodes_have_descendant_nodes(self):
        self.assertEqual(set(self.a1.descendant_nodes),
                         set([self.a3, self.a4, self.a5, self.a6]))

    def test_descendant_nodes_have_ancestor_nodes(self):
        self.assertEqual(set(self.a6.ancestor_nodes), set([self.a1, self.a2]))


class QueryParentChildTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A(a_name='a1')
        self.b1 = B(b_name='b1', parent_nodes=[self.a1])
        self.a2 = A(a_name='a2', parent_nodes=[self.b1])
        self.b2 = B(b_name='b2', parent_nodes=[self.a2])
        self.session.add_all([self.a1, self.b1, self.a2, self.b2])

    def test_filter_by_parent_attrs(self):
        q = (
            self.session.query(B)
            .join(B.parent_nodes, aliased=True)
            .join(A)
            .filter(A.a_name == 'a1')
        )
        self.assertEqual(q.all(), [self.b1])

    def test_filter_by_child_attrs(self):
        q = (
            self.session.query(A)
            .join(A.child_nodes, aliased=True)
            .join(B)
            .filter(B.b_name == 'b2')
        )
        self.assertEqual(q.all(), [self.a2])


class QueryScionsAndParentsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A(a_name='a1')
        self.b1 = B(b_name='b1', parent_nodes=[self.a1],
                    ancestor_nodes=[self.a1])
        self.a2 = A(a_name='a2', parent_nodes=[self.b1],
                    ancestor_nodes=[self.a1])
        self.b2 = B(b_name='b2', parent_nodes=[self.a2],
                    ancestor_nodes=[self.a1])
        self.b3 = B(b_name='b3', parent_nodes=[self.a1],
                    ancestor_nodes=[self.a1])
        self.session.add_all([self.a1, self.b1, self.a2, self.b2])

    def test_can_filter_by_parent_and_ancestor(self):
        q = (
            self.session.query(B)
            .join(B.ancestor_nodes, aliased=True)
            .join(A, aliased=True)
            .filter(A.a_name == 'a1')
            .join(B.parent_nodes, aliased=True, from_joinpoint=True)
            .join(A, aliased=True)
            .filter(A.a_name == 'a2')
        )
        self.assertEqual(q.all(), [self.b2])


class QueryChildrenParentsOfTypeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.a1 = A()
        self.b1 = B()
        self.a2 = A(parent_nodes=[self.a1, self.b1])
        self.b2 = B(parent_nodes=[self.a1, self.b1])
        self.session.add_all([self.a1, self.b1, self.a2, self.b2])

    def test_child_nodes_of_type(self):
        self.assertEqual(self.a1.child_nodes_of_type('A'), [self.a2])
        self.assertEqual(self.a1.child_nodes_of_type('B'), [self.b2])

    def test_parent_nodes_of_type(self):
        self.assertEqual(self.a2.parent_nodes_of_type('A'), [self.a1])
        self.assertEqual(self.a2.parent_nodes_of_type('B'), [self.b1])
