import collections
import importlib


class DotSpecLoader(object):
    @classmethod
    def load_from_dot_spec(self, dot_spec=None):
        dot_spec_parts = dot_spec.split(':')
        if len(dot_spec_parts) > 2: raise Exception("invalid dot spec")
        module_name = dot_spec_parts[0]
        module = importlib.import_module(module_name)
        result = module
        if len(dot_spec_parts) == 2:
            object_dot_spec = dot_spec_parts[1]
            result = self.get_obj_value_from_dot_spec(obj=module,
                                                      dot_spec=object_dot_spec)
        return result

    @classmethod
    def get_obj_value_from_dot_spec(self, obj=None, dot_spec=None):
        path_elements = dot_spec.split('.')
        cursor = obj
        pos = 0
        try:
            for path_element in path_elements:
                next_cursor = self.get_attr_or_item(obj=cursor, key=path_element)
                cursor = next_cursor
                pos += 1
            return cursor
        except Exception as exception:
            raise Exception(("Unable to get value\n"
                             "  exception: {exception}\n"
                             "  error at path_element {pos} ({partial_path})"
                            ).format(exception=exception,
                                     pos=pos,
                                     partial_path=path_elements[:pos+1]))

    @classmethod
    def get_attr_or_item(self, obj=None, key=None):
        try:
            if isinstance(obj, collections.abc.Sequence): return obj[int(key)]
            elif isinstance(obj, collections.abc.Mapping): return obj[key]
            elif hasattr(obj, key): return getattr(obj, key)
        except: return None
