import unittest

from mc.db.db import Db
from .. import basic_entity_selector


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.db = Db(db_uri='sqlite://', ensure_tables=True)
        self.tags = ['a', 'b', 'c']
        self.ent_type = 'some_ent_type'
        self.objects = self.setup_objects()
        self.selector = basic_entity_selector.BasicEntitySelector(db=self.db)

    def setup_objects(self):
        ents = self.setup_ents()
        return {
            'ents': ents,
            'requests': self.setup_requests(ents=ents)
        }

    def setup_ents(self):
        tag_combos = [[], ['a'], ['a', 'b'], ['a', 'b', 'c']]
        ents = [
            self.db.models.Ent(
                ent_type=self.ent_type,
                tags=tag_combo
            )
            for tag_combo in tag_combos
        ]
        self.db.session.add_all(ents)
        self.db.session.commit()
        return ents

    def setup_requests(self, ents=None):
        requests = []
        requests.extend([
            self.db.models.Request(
                request_type='some_request_type',
                request_tag='request_tag_1',
                instance_key=ent.key
            )
            for ent in [ents[1]]
        ])
        requests.extend([
            self.db.models.Request(
                request_type='some_request_type',
                request_tag='request_tag_2',
                instance_key=ent.key
            )
            for ent in [ents[2], ents[3]]
        ])
        requests.append(self.db.models.Request(
            request_type='some_request_type',
            request_tag='lacks_parents',
            instance_key='some_instance_key'
        ))
        self.db.session.add_all(requests)
        self.db.session.commit()
        return requests

    def ents_to_expected_items(self, ents=None):
        return [
            {'key': ent.key, 'value': ent}
            for ent in sorted(ents, key=(lambda _ent: _ent.modified))
        ]


class UnfilteredTestCase(BaseTestCase):
    def test_gets_all_items(self):
        items = list(self.selector.get_items(ent_type=self.ent_type))
        expected_items = self.ents_to_expected_items(ents=self.objects['ents'])
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class HavingTagsTestCase(BaseTestCase):
    def test_gets_items_w_tags(self):
        tags = self.tags[:2]
        items = list(self.selector.get_items(ent_type=self.ent_type,
                                             having_tags=tags))
        expected_items = self.ents_to_expected_items(
            ents=[ent for ent in self.objects['ents']
                  if set(ent.tags) >= set(tags)]
        )
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class LackingTagsTestCase(BaseTestCase):
    def test_gets_items_sans_tags(self):
        tags_to_exclude = self.tags[2:]
        items = list(self.selector.get_items(
            ent_type=self.ent_type, sans_tags=tags_to_exclude))
        expected_items = self.ents_to_expected_items(
            ents=[ent for ent in self.objects['ents']
                  if not (set(ent.tags) & set(tags_to_exclude))]
        )
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class LackingRequestsWTagsTestCase(BaseTestCase):
    def test_gets_items_sans_matching_requests(self):
        tags_to_exclude = ['request_tag_1']
        items = list(self.selector.get_items(
            ent_type=self.ent_type,
            sans_requests_w_tags=tags_to_exclude)
        )
        expected_items = self.ents_to_expected_items(
            ents=[
                ent for ent in self.objects['ents']
                if (
                    not (
                        set(tags_to_exclude) &
                        set([request.request_tag
                             for request in self.objects['requests']
                             if request.instance_key == ent.key])
                    )
                )
            ]
        )
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class CombinedFiltersTestCase(BaseTestCase):
    def test_gets_items_matching_combined_filters(self):
        having_tags = ['a', 'b']
        sans_tags = ['c']
        sans_requests_w_tags = ['request_tag_1']
        items = list(self.selector.get_items(
            ent_type=self.ent_type,
            having_tags=having_tags,
            sans_tags=sans_tags,
            sans_requests_w_tags=sans_requests_w_tags
        ))

        filters = []

        def having_tags_filter(ent): return set(ent.tags) >= set(having_tags)
        filters.append(having_tags_filter)

        def sans_tags_filter(ent):
            return not (set(ent.tags) & set(sans_tags))
        filters.append(sans_tags_filter)

        def sans_requests_w_tags_filter(ent):
            related_request_tags = [
                req.request_tag for req in self.objects['requests']
                if (req.instance_key == ent.key)
            ]
            return not (
                set(related_request_tags) & set(sans_requests_w_tags)
            )
        filters.append(sans_requests_w_tags_filter)

        expected_ents = [ent for ent in self.objects['ents']
                         if all(filter_(ent) for filter_ in filters)]
        expected_items = self.ents_to_expected_items(ents=expected_ents)
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)
