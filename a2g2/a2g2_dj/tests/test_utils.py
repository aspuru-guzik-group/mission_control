from unittest.mock import call, MagicMock, patch

from django.test import TestCase

from .. import utils as _utils


@patch.object(_utils, 'ChemThingActionProcessor')
@patch.object(_utils, 'deserialize_bulk_actions')
class ProcessSerializedChemThingActionsTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.serialized_actions = MagicMock()

    def _do_process(self):
        return _utils.process_serialized_chemthing_actions(
            serialized_chemthing_actions=self.serialized_actions)

    def test_deserializes(self, *args):
        self._do_process()
        self.assertEqual(_utils.deserialize_bulk_actions.call_args,
                         call(serialized_bulk_actions=self.serialized_actions))

    def test_passes_deserialized_actions_action_processor(self, *args):
        result = self._do_process()
        expected_processor = _utils.ChemThingActionProcessor.return_value
        self.assertEqual(
            expected_processor.process_actions.call_args,
            call(actions=_utils.deserialize_bulk_actions.return_value))
        self.assertEqual(result,
                         expected_processor.process_actions.return_value)
