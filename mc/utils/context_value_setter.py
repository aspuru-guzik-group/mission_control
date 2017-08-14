"""
Utils for gettings/setting values via python dot paths.
Also includes logic for doing transforms on values during get/set operations.
"""

import collections
import copy
import json

from . import dot_spec_loader


class ContextValueSetter(object):
    class UnknownTransformError(Exception): pass
    class InvalidValueSpecError(Exception): pass
    class InvalidMappingSourceError(Exception):
        def __init__(self, *args, mapping_source=None, **kwargs):
            error = "invalid mapping_source '{}'".format(mapping_source)
            super().__init__(error, *args, **kwargs)

    def set_context_values(self, value_specs=None, context=None):
        for value_spec in (value_specs or []):
            self.set_context_value(value_spec=value_spec, context=context)

    def set_context_value(self, value_spec=None, context=None):
        """
        Args:
            value_spec (dict, optional): a dict with this shape:
                ::

                    {
                        value|source: a value or source string,
                        transform (optional): an optional transform to apply.
                            See :module:`transform_value_for_value_spec`
                        dest: a dot string specifying where to put
                            the set value
                    }

            context (dict): an object to use as the 'ctx' root for dot strings.
        """
        try:
            value_spec = self.compile_value_spec(value_spec=value_spec)
            value = self.get_value_for_value_spec(
                value_spec=value_spec, context=context)
            transformed_value = self.transform_value_for_value_spec(
                value=value, value_spec=value_spec, context=context)
            self.set_context_dest_value(
                context=context, dest=value_spec['dest'], 
                value=transformed_value)
        except Exception as exception:
            raise Exception("Unable to set_context_value:\n"
                            "  exception: {exception}\n"
                            "  value_spec: '{value_spec}'\n"
                            "  context: '{context}'".format(
                                exception=exception,
                                value_spec=value_spec,
                                context=self.elide_string(str(context))
                            ))

    def compile_value_spec(self, value_spec):
        if type(value_spec) is str:
            try:
                source, dest = value_spec.split('=>')
                value_spec = {'source': source.strip(), 'dest': dest.strip()}
            except: raise self.InvalidValueSpecError(value_spec)
        return value_spec

    def elide_string(self, string=None, max_len=100):
        elided_string = string
        if len(string) > max_len:
            elided_string = string[:max_len] + '...'
        return elided_string

    def get_value_for_value_spec(self, value_spec=None, context=None):
        if 'value' in value_spec: value = value_spec['value']
        elif 'source' in value_spec:
            value = self.get_value_from_dot_spec(obj={'ctx': context},
                                                 dot_spec=value_spec['source'])
        else: 
            raise Exception("Can't handle value_spec '{value_spec}'".format(
                value_spec=value_spec))
        return value

    def transform_value_for_value_spec(self, value=None, value_spec=None,
                                       context=None):
        """Transform a value per a value spec's 'transform' prop.

        Transform types can be:
            - json.dumps
            - json.loads
            - mapping (see :meth:`execute_mapping_transform`)
        """
        transform = value_spec.get('transform')
        if not transform: return value
        if type(transform) is str: transform = {'type': transform}
        if transform['type'] == 'json.dumps': return json.dumps(value)
        elif transform['type'] == 'json.loads': return json.loads(value)
        elif transform['type'] == 'mapping':
            return self.execute_mapping_transform(
                value=value, transform=transform, context=context)
        raise self.UnknownTransformError(transform)

    def execute_mapping_transform(self, value=None, transform=None,
                                  context=None):
        """
        Args:
            value (obj): value to use items source for mapping.
            transform (dict, optional): transform spec, which is a dict
                like this: ::

                    {
                        type: 'transform',
                        params: {
                            <skeleton>: <a template for the mapping product>
                            [wirings]: a list of wirings to perform on the
                            skeleton.
                        }
                    }

                The  wirings will be the given context plus an 'item' prop
                representing the item being mapped. This prop will be a dict
                with keys 'item', 'items', and 'idx'

        Returns:
            mapped_items (list): the mapped items
        """
        iterator = self.get_iterator_for_mapping_source(mapping_source=value)
        mapped_items = [
            self.map_item(
                item={'key': item_key, 'value': item_value, 'idx': item_idx},
                items=value, mapping_params=transform['params'],
                context=context
            ) for item_idx, (item_key, item_value) in enumerate(iterator)
        ]
        return mapped_items

    def get_iterator_for_mapping_source(self, mapping_source=None):
        iterator = None
        if isinstance(mapping_source, collections.abc.Mapping):
            iterator = mapping_source.items()
        elif isinstance(mapping_source, collections.abc.Sequence):
            iterator = enumerate(mapping_source)
        if not iterator: raise self.InvalidMappingSourceError(mapping_source)
        return iterator

    def map_item(self, item=None, items=None, mapping_params=None, context=None):
        skeleton_copy = copy.deepcopy(mapping_params.get('skeleton', {}))
        self.set_context_values(
            value_specs=mapping_params.get('wirings', []),
            context={**context, 'item': item,
                     'items': items, 'skeleton': skeleton_copy})
        return skeleton_copy

    def set_context_dest_value(self, context=None, dest=None, value=None):
        self.set_value_from_dot_spec(obj={'ctx': context}, dot_spec=dest,
                                     value=value)

    def get_attr_or_item(self, obj=None, key=None):
        return dot_spec_loader.DotSpecLoader.get_attr_or_item(obj=obj, key=key)

    def set_value_from_dot_spec(self, obj=None, value=None, dot_spec=None):
        path_elements = dot_spec.split('.')
        cursor = obj
        pos = 0
        try:
            for path_element in path_elements[:-1]:
                next_cursor = self.get_attr_or_item(obj=cursor,
                                                    key=path_element)
                if next_cursor is None:
                    cursor[path_element] = {}
                    next_cursor = cursor[path_element]
                cursor = next_cursor
                pos += 1
            last_key = path_elements[-1]
            if self.is_listlike(cursor): last_key = int(last_key)
            cursor[last_key] = value

        except Exception as exception:
            raise Exception(("Unable to set value\n"
                             "  exception: {exception}\n"
                             "  error at path_element {pos} ({partial_path})"
                            ).format(exception=exception,
                                     pos=pos,
                                     partial_path=path_elements[:pos+1]))

    def get_value_from_dot_spec(self, obj=None, dot_spec=None):
        return dot_spec_loader.DotSpecLoader.get_obj_value_from_dot_spec(
            obj=obj, dot_spec=dot_spec)

    def is_listlike(self, o):
        return isinstance(o, collections.Sequence) and not isinstance(o, str)


def set_context_values(value_specs=None, context=None):
    ContextValueSetter().set_context_values(value_specs=value_specs,
                                            context=context)

def set_context_value(value_spec=None, context=None):
    ContextValueSetter().set_context_value(value_spec=value_spec,
                                           context=context)

def execute_mapping_transform(items=None, transform=None, context=None):
    return ContextValueSetter().execute_mapping_transform(
        value=items, transform=transform, context=context)
