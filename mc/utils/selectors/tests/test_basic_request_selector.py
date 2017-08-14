import unittest

from mc.db.db import Db
from .. import basic_request_selector


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.db = Db(db_uri='sqlite://', ensure_tables=True)
        self.selector = basic_request_selector.BasicRequestSelector(db=self.db)
        self.reqs = self.setup_reqs()

    def setup_reqs(self):
        reqs = [
            self.db.models.Request(
                status='status_%s' % i,
                request_tag='req_tag_%s' % i,
                request_type='req_type_%s' % i,
                instance_key=str(i)
            )
            for i in range(3)
        ]
        self.db.session.add_all(reqs)
        self.db.session.commit()
        return reqs

    def reqs_to_expected_items(self, reqs=None):
        return [
            {'key': req.key, 'value': req}
            for req in sorted(reqs, key=(lambda _req: _req.modified))
        ]


class UnfilteredTestCase(BaseTestCase):
    def test_gets_all_requests(self):
        items = list(self.selector.get_items())
        expected_items = self.reqs_to_expected_items(reqs=self.reqs)
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class HavingRequestTagTestCase(BaseTestCase):
    def test_gets_items_w_request_tag(self):
        items = list(self.selector.get_items(
            having_request_tag=self.reqs[0].request_tag))
        expected_items = self.reqs_to_expected_items(reqs=[self.reqs[0]])
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class HavingRequestTypeTestCase(BaseTestCase):
    def test_gets_items_w_request_type(self):
        items = list(self.selector.get_items(
            having_request_type=self.reqs[0].request_type))
        expected_items = self.reqs_to_expected_items(reqs=[self.reqs[0]])
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class HavingStatusTestCase(BaseTestCase):
    def test_gets_items_w_status(self):
        items = list(self.selector.get_items(
            having_status=self.reqs[0].status))
        expected_items = self.reqs_to_expected_items(reqs=[self.reqs[0]])
        self.assertTrue(len(items) > 0)
        self.assertEqual(items, expected_items)


class LackingDerivedRequestsWTagsTestCase(BaseTestCase):
    def setup_reqs(self):
        reqs = [
            self.db.models.Request(
                request_type='some_request_type',
                request_tag='req_tag_%s' % i,
                instance_key=str(i)
            )
            for i in range(3)
        ]
        self.db.session.add_all(reqs)
        self.db.session.commit()
        child_reqs = [
            self.db.models.Request(
                request_type='some_request_type',
                instance_key=parent_req.key,
                request_tag=(parent_req.request_tag + '_child')
            )
            for parent_req in reqs
        ]
        self.db.session.add_all(child_reqs)
        self.db.session.commit()
        return reqs + child_reqs

    def test_gets_items_sans_downstream_requests(self):
        tags_to_exclude = ['req_tag_1_child']
        items = list(self.selector.get_items(
            sans_downstream_requests_w_tags=tags_to_exclude))
        expected_items = self.reqs_to_expected_items(
            reqs=[
                req for req in self.reqs
                if (
                    not (
                        set(tags_to_exclude) &
                        set([req_.request_tag for req_ in self.reqs
                             if req_.instance_key == req.key])
                    )
                )
            ]
        )
        self.assertTrue(len(items) > 0)
        self.assertTrue(len(items) < len(self.reqs))
        self.assertEqual(items, expected_items)
