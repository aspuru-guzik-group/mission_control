from unittest.mock import call, MagicMock, patch

from django.test import TestCase

from .. import chemthing_action_processor


class BaseTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.processor = chemthing_action_processor.ChemThingActionProcessor()

class ProcessActionsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.actions = list(range(3))
        self.processor.process_action = MagicMock()
        self.processor.process_actions(actions=self.actions)

    def test_calls_process_action_for_each_item(self):
        expected_call_args_list = [call(action=action)
                                   for action in self.actions]
        self.assertEqual(self.processor.process_action.call_args_list,
                         expected_call_args_list)

class ProcessActionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.action = MagicMock()
        for attr in ['get_or_create_chemthing_for_action',
                     'update_chemthing_from_action']:
            setattr(self.processor, attr, MagicMock())
        self.result = self.processor.process_action(action=self.action)

    def test_updates_get_or_create_result(self):
        self.assertEqual(
            self.processor.get_or_create_chemthing_for_action.call_args,
            call(action=self.action))
        expected_get_or_create_result = \
                self.processor.get_or_create_chemthing_for_action.return_value
        self.assertEqual(
            self.processor.update_chemthing_from_action.call_args,
            call(chemthing=expected_get_or_create_result, action=self.action)) 
        self.assertEqual(
            self.result,
            self.processor.update_chemthing_from_action.return_value)

class GetOrCreateChemForActionThing(BaseTestCase):
    @patch.object(chemthing_action_processor, 'ChemThing')
    def test_dispatches_to_get_or_create(self, MockChemThing):
        action_with_key = {'key': MagicMock()}
        self.processor.get_or_create_chemthing_for_action(
            action=action_with_key)
        self.assertEqual(MockChemThing.objects.get_or_create.call_args,
                         call(keys__has_key=action_with_key['key']))
        action_sans_key = {}
        self.processor.get_or_create_chemthing_for_action(
            action=action_sans_key)
        self.assertEqual(MockChemThing.objects.get_or_create.call_args, call())

class UpdateChemThingFromAction(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.action = {
            'updates': {attr: MagicMock() for attr in (
                self.processor.ATTRS_TO_UPDATE + self.processor.ATTRS_TO_SET)}
        }
        self.chemthing = MagicMock()
        self.processor.update_chemthing_from_action(chemthing=self.chemthing,
                                                    action=self.action)

    def test_processes_attrs_to_update(self):
        expected_values = {}
        actual_values = {}
        for attr in self.processor.ATTRS_TO_UPDATE:
            actual_values[attr] = getattr(self.chemthing, attr)
            expected_values[attr] = {**getattr(self.chemthing, attr, {}),
                                     **self.action['updates'][attr]}
        self.assertEqual(actual_values, expected_values)

    def test_processes_attrs_to_set(self):
        expected_values = {}
        actual_values = {}
        for attr in self.processor.ATTRS_TO_SET:
            actual_values[attr] = getattr(self.chemthing, attr)
            expected_values[attr] = self.action['updates'][attr]
        self.assertEqual(actual_values, expected_values)

    def test_saves_chemthing(self):
        self.assertEqual(self.chemthing.save.call_args, call())
