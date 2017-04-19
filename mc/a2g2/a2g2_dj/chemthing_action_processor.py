from .models import ChemThing


class ChemThingActionProcessor(object):
    ATTRS_TO_UPDATE = ['props', 'types', 'keys', 'precursors', 'ancestors']
    ATTRS_TO_SET = ['uuid']

    def process_actions(self, actions=None):
        for action in actions: self.process_action(action=action)

    def process_action(self, action=None):
        chemthing = self.get_or_create_chemthing_for_action(action=action)
        return self.update_chemthing_from_action(chemthing=chemthing, action=action)

    def get_or_create_chemthing_for_action(self, action=None):
        filter_kwargs = {}
        if 'key' in action: filter_kwargs['key'] = action['key']
        return ChemThing.objects.get_or_create(**filter_kwargs)

    def update_chemthing_from_action(self, chemthing=None, action=None):
        updates = action.get('updates', {})
        for attr in self.ATTRS_TO_UPDATE:
            if attr in updates:
                prev_value = getattr(chemthing, attr, {})
                next_value = {**prev_value, **updates[attr]}
                setattr(chemthing, attr, next_value)
        for attr in self.ATTRS_TO_SET:
            if attr in updates: setattr(chemthing, attr, updates[attr])
        chemthing.save()
