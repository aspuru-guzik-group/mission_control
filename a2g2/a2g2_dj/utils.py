import json

from .chemthing_action_processor import ChemThingActionProcessor


def process_serialized_chemthing_actions(serialized_chemthing_actions=None):
    processor = ChemThingActionProcessor()
    actions = deserialize_bulk_actions(
        serialized_bulk_actions=serialized_chemthing_actions)
    return processor.process_actions(actions=actions)

def deserialize_bulk_actions(serialized_bulk_actions=None):
    for serialized_action in serialized_bulk_actions.split("\n"):
        yield json.loads(serialized_action)
