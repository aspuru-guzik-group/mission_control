import time
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
        self.ent_type = 'some_ent_type'
        self.props = {
            'str_prop': 'str',
            'int_prop': 1,
            'bool_prop': True,
            'mapping_prop': {'some': {'nested': 'prop'}},
            'sequence_prop': ['some', 'sequence']
        }
        self.tags = {'tag_%s' % i for i in range(3)}

    def generate_ent(self, **kwargs):
        return models.Ent(**{
            'ent_type': self.ent_type,
            'props': self.props,
            'tags': self.tags,
            **kwargs
        })


class CreateEntTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.ent = self.generate_ent()
        self.session.add(self.ent)
        self.session.commit()

    def test_has_timestamps(self):
        self.assertTrue(self.ent.created is not None)
        self.assertTrue(self.ent.modified is not None)

    def test_default_key_startswith_ent_ent_type(self):
        self.assertTrue(self.ent.key.startswith('ent:' + self.ent_type))

    def test_has_props(self):
        self.assertTrue(self.ent.props, self.props)

    def test_has_tags(self):
        self.assertTrue(self.ent.tags, self.tags)


class QueryEntTestCase(BaseTestCase):
    def test_filter_created(self):
        ents = []
        for i in range(3):
            ent = self.generate_ent()
            self.session.add(ent)
            self.session.commit()
            ents.append(ent)
            time.sleep(1e-3)
        newer_than_ent_0 = (self.session.query(models.Ent)
                            .filter(models.Ent.created > ents[0].created)
                            .all())
        self.assertEqual(newer_than_ent_0, ents[1:])
        older_than_ent_2 = (self.session.query(models.Ent)
                            .filter(models.Ent.created < ents[2].created)
                            .all())
        self.assertEqual(older_than_ent_2, ents[:2])

    def test_filter_props(self):
        ents = []
        for i in range(3):
            ent = self.generate_ent(
                props={'int_prop': i, 'str_prop': 'str_%s' % i}
            )
            self.session.add(ent)
            self.session.commit()
            ents.append(ent)
        ents_w_int_eq_1 = (
            self.session.query(models.Ent)
            .filter(models.Ent.props_set.any(key='int_prop', value=1))
            .all()
        )
        self.assertEqual(ents_w_int_eq_1, [ents[1]])
        ents_w_int_gt_1 = (
            self.session.query(models.Ent)
            .filter(
                models.Ent.props_set.any(
                    (models.Ent.Prop.key == 'int_prop')
                    & (models.Ent.Prop.value > 1)
                )
            )
            .all()
        )
        self.assertEqual(ents_w_int_gt_1, ents[2:])

    def test_filter_tags(self):
        ents = [
            self.generate_ent(tags={'tag_1', 'tag_2'}),
            self.generate_ent(tags={'tag_2', 'tag_3'}),
            self.generate_ent(tags={'tag_1', 'tag_3'}),
        ]
        self.session.add_all(ents)
        self.session.commit()
        ents_w_tag_1 = (
            self.session.query(models.Ent)
            .filter(models.Ent.tags_set.any(name='tag_1'))
            .all()
        )
        self.assertEqual(set(ents_w_tag_1), set([ents[0], ents[2]]))
        ents_w_tag_2 = (
            self.session.query(models.Ent)
            .filter(models.Ent.tags_set.any(name='tag_2'))
            .all()
        )
        self.assertEqual(set(ents_w_tag_2), set([ents[0], ents[1]]))
        ents_w_nonexistent_tag = (
            self.session.query(models.Ent)
            .filter(models.Ent.tags_set.any(name='nonexistent_tag'))
            .all()
        )
        self.assertEqual(set(ents_w_nonexistent_tag), set())


class LineageTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.families = self._create_families()

    def _create_families(self):
        families = {}
        for i in range(2):
            family_key = ('family_%s' % i)
            families[family_key] = self._create_family(family_key=family_key)
            return families

    def _create_family(self, family_key=None):
        common_props = {'family_key': family_key}
        grandparents = [
            self.generate_ent(
                key=('%s:grandparent_%s' % (family_key, i)),
                props={**common_props, 'generation': 'grandparents'}
            )
            for i in range(4)
        ]
        grandparent_pairs = [
            [grandparents[0], grandparents[1]],
            [grandparents[2], grandparents[3]]
        ]
        parents = []
        for i, grandparent_pair in enumerate(grandparent_pairs):
            parents.append(
                self.generate_ent(
                    key=('%s:parent_%s' % (family_key, i)),
                    props={**common_props, 'generation': 'parents'},
                    parents=grandparent_pair,
                    ancestors=grandparent_pair
                )
            )
        children = [
            self.generate_ent(
                key=('%s:child_%s' % (family_key, i)),
                props={**common_props, 'generation': 'children'},
                parents=parents,
                ancestors=(grandparents + parents)
            )
            for i in range(3)
        ]
        self.session.add_all(grandparents + parents + children)
        self.session.commit()
        family = {
            'grandparents': grandparents,
            'grandparent_pairs': grandparent_pairs,
            'parents': parents,
            'children': children
        }
        return family

    def test_parents(self):
        for family in self.families.values():
            for child in family['children']:
                self.assertEqual(child.parents, family['parents'])
                for i, parent in enumerate(family['parents']):
                    self.assertEqual(
                        parent.parents,
                        family['grandparent_pairs'][i]
                    )

    def test_children(self):
        for family in self.families.values():
            for parent in family['parents']:
                self.assertEqual(
                    set(parent.children),
                    set(family['children'])
                )
                for i, gp_pair in enumerate(family['grandparent_pairs']):
                    for grandparent in gp_pair:
                        self.assertEqual(
                            set(grandparent.children),
                            set([family['parents'][i]])
                        )

    def test_ancestors(self):
        for family in self.families.values():
            for child in family['children']:
                self.assertEqual(
                    set(child.ancestors),
                    set(family['grandparents'] + family['parents'])
                )

    def test_descendants(self):
        for family in self.families.values():
            for grandparent in family['grandparents']:
                self.assertEqual(
                    set(grandparent.descendants),
                    set(grandparent.children + family['children'])
                )

    def test_query_on_parents(self):
        for family in self.families.values():
            children_of_grandparent_pair_0 = (
                self.session.query(models.Ent)
                .join(models.Ent.parents, aliased=True, from_joinpoint=True)
                .filter(
                    models.Ent.key.in_([
                        grandparent.key
                        for grandparent in family['grandparent_pairs'][0]
                    ])
                )
                .reset_joinpoint()
                .all()
            )
            self.assertEqual(children_of_grandparent_pair_0,
                             [family['parents'][0]])

    def test_query_on_ancestors(self):
        for family_key, family in self.families.items():
            descendants = (
                self.session.query(models.Ent)
                .filter(
                    models.Ent.props_set.any(
                        key='generation',
                        value='children'
                    )
                )
                .join(models.Ent.ancestors, aliased=True, from_joinpoint=True)
                .filter(
                    models.Ent.props_set.any(
                        key='family_key', value=family_key
                    )
                )
                .reset_joinpoint()
                .all()
            )
            self.assertEqual(set(descendants), set(family['children']))
