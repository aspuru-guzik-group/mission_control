import sqlalchemy as _sqla

from mc.utils import hash_utils
from . import utils


common_supers = [utils.KeyedMixin, utils.TimestampMixin, utils.PropsMixin,
                 utils.TagsMixin, utils.ToDictMixin, utils.Base]


class Flow(*common_supers):
    claimed = utils.generate_boolean_column()
    cfg = utils.generate_json_column()
    data = utils.generate_json_column()
    graph = utils.generate_json_column()
    num_tickable_tasks = utils.generate_int_column()
    depth = utils.generate_int_column()


class Job(*common_supers):
    HASH_COMPONENTS = ['params', 'job_type']
    claimed = utils.generate_boolean_column()
    status = utils.generate_status_column()
    job_type = utils.generate_str_column()
    params = utils.generate_json_column()
    cfg = utils.generate_json_column()
    data = utils.generate_json_column()
    artifact_meta = utils.generate_json_column()
    hash = utils.generate_str_column(length=63, unique=True, nullable=True)

    @classmethod
    def _receive_hash_component(cls, component_name=None, target=None,
                                value=None, oldvalue=None, **kwargs):
        if value == oldvalue:
            return
        if component_name == 'module' and value is None:
            return
        components = {
            component_name: value,
            **{
                other_component_name: getattr(target, other_component_name)
                for other_component_name in Job.HASH_COMPONENTS
                if other_component_name != component_name
            }
        }
        target.hash = hash_utils.hash_obj(components)


for component_name in Job.HASH_COMPONENTS:
    def receiver(component_name=component_name, **kwargs):
        return Job._receive_hash_component(component_name=component_name,
                                           **kwargs)
    listener_decorator = _sqla.event.listens_for(
        getattr(Job, component_name), 'set', named=True)
    listener_decorator(receiver)


class Queue(*common_supers):
    queue_spec = utils.generate_json_column()


class Lock(*common_supers):
    lockee_key = utils.generate_key_column('lockee_key', default=None)
    locker_key = utils.generate_key_column('locker_key', default=None)
