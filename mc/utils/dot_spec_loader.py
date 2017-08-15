"""Utilities for loading values from python dot specs paths.
"""

import collections
import importlib


class DotSpecLoader(object):
    @classmethod
    def load_from_dot_spec(cls, dot_spec=None):
        """
        Args:
            dot_spec (str): dot_spec should be a period-separated path,
                like a python module path.  A dot_spec can two parts separated
                by a ':'. In this case the first part of the path will be
                interpreted as a module path, and the second part of the path
                will be interpreted as the dot_path to a member of that module.
        """
        return cls()._load_from_dot_spec(dot_spec=dot_spec)

    def _load_from_dot_spec(self, dot_spec=None):
        dot_spec_parts = dot_spec.split(':')
        if len(dot_spec_parts) > 2:
            raise Exception("invalid dot spec")
        module_name = dot_spec_parts[0]
        module = importlib.import_module(module_name)
        result = module
        if len(dot_spec_parts) == 2:
            object_dot_spec = dot_spec_parts[1]
            result = self._get_obj_value_from_dot_spec(
                obj=module, dot_spec=object_dot_spec)
        return result

    @classmethod
    def get_obj_value_from_dot_spec(cls, obj=None, dot_spec=None):
        return cls()._get_obj_value_from_dot_spec(obj=obj, dot_spec=dot_spec)

    def _get_obj_value_from_dot_spec(self, obj=None, dot_spec=None):
        path_elements = dot_spec.split('.')
        cursor = obj
        pos = 0
        try:
            for path_element in path_elements:
                next_cursor = self._get_attr_or_item(obj=cursor,
                                                     key=path_element)
                cursor = next_cursor
                pos += 1
            return cursor
        except Exception as exception:
            raise Exception(
                (
                    "Unable to get value\n"
                    "  exception: {exception}\n"
                    "  error at path_element {pos} ({partial_path})"
                ).format(exception=exception, pos=pos,
                         partial_path=path_elements[:pos+1])
            )

    @classmethod
    def get_attr_or_item(cls, obj=None, key=None):
        return cls()._get_attr_or_item(obj=obj, key=key)

    def _get_attr_or_item(self, obj=None, key=None):
        if isinstance(key, str) and key.endswith('()'):
            return self._get_attr_or_item(obj, key[:-2])()
        for type_ in ['sequence', 'mapping', 'mapping_view', 'attr']:
            try:
                return getattr(self, '_get_from_' + type_)(obj, key)
            except:
                pass
        return None

    def _get_from_sequence(self, obj=None, key=None):
        if isinstance(obj, collections.abc.Sequence):
            return obj[int(key)]
        raise Exception("Not a sequence")

    def _get_from_mapping(self, obj=None, key=None):
        if isinstance(obj, collections.abc.Mapping):
            return obj[key]
        raise Exception("Not a mapping")

    def _get_from_mapping_view(self, obj=None, key=None):
        if isinstance(obj, collections.abc.MappingView):
            return list(obj)[int(key)]
        raise Exception("Not a mapping_view")

    def _get_from_attr(self, obj=None, key=None):
        if hasattr(obj, key):
            return getattr(obj, key)
        raise Exception("Does not have attr")


def load_from_dot_spec(dot_spec=None):
    return DotSpecLoader.load_from_dot_spec(dot_spec=dot_spec)


def get_attr_or_item(obj=None, key=None):
    return DotSpecLoader.get_attr_or_item(obj=obj, key=key)
