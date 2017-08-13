import sqlalchemy.orm as _orm
from mc.utils import hash_utils

from . import constants
from . import utils


class Flow(*utils.common_supers):
    label = utils.generate_str_column()
    claimed = utils.generate_boolean_column()
    status = utils.generate_status_column()
    parent_key = utils.generate_str_column(length=constants.KEY_LENGTH)
    cfg = utils.generate_json_column()
    data = utils.generate_json_column()
    graph = utils.generate_json_column()
    num_tickable_tasks = utils.generate_int_column()
    depth = utils.generate_int_column()
    __mapper_args__ = {
        'polymorphic_identity': 'flow',
    }


class Job(*utils.common_supers):
    HASH_COMPONENTS = ['job_type', 'job_params']
    label = utils.generate_str_column()
    claimed = utils.generate_boolean_column()
    status = utils.generate_status_column()
    parent_key = utils.generate_str_column(length=constants.KEY_LENGTH)
    job_type = utils.generate_str_column()
    job_params = utils.generate_json_column()
    cfg = utils.generate_json_column()
    data = utils.generate_json_column()
    artifact_meta = utils.generate_json_column()
    job_hash = utils.generate_str_column(length=63, unique=True,
                                         nullable=True)
    __mapper_args__ = {
        'polymorphic_identity': 'job',
    }

    @classmethod
    def _receive_hash_component(cls, component_name=None, target=None,
                                value=None, oldvalue=None, **kwargs):
        if value == oldvalue:
            return
        if component_name == 'job_type' and value is None:
            return
        components = {
            component_name: value,
            **{
                other_component_name: getattr(target, other_component_name)
                for other_component_name in Job.HASH_COMPONENTS
                if other_component_name != component_name
            }
        }
        target.job_hash = hash_utils.hash_obj(components)


for component_name in Job.HASH_COMPONENTS:
    def receiver(component_name=component_name, **kwargs):
        return Job._receive_hash_component(component_name=component_name,
                                           **kwargs)
    listener_decorator = utils._sqla.event.listens_for(
        getattr(Job, component_name), 'set', named=True)
    listener_decorator(receiver)


class Queue(*utils.common_supers):
    label = utils.generate_str_column()
    queue_spec = utils.generate_json_column()
    __mapper_args__ = {
        'polymorphic_identity': 'queue',
    }


class Lock(*utils.common_supers):
    lockee_key = utils.generate_str_column(length=constants.KEY_LENGTH)
    locker_key = utils.generate_str_column(length=constants.KEY_LENGTH)
    __mapper_args__ = {
        'polymorphic_identity': 'lock',
    }


class Request(*utils.common_supers):
    HASH_COMPONENTS = ['request_type', 'instance_key', 'params']
    request_type = utils.generate_str_column(primary_key=True)
    request_tag = utils.generate_str_column(primary_key=True)
    instance_key = utils.generate_str_column(primary_key=True)
    params = utils.generate_json_column()
    status = utils.generate_status_column(default='PENDING')
    __mapper_args__ = {
        'polymorphic_identity': 'request',
    }


class Ent(*utils.common_supers):
    ent_type = utils.generate_str_column(
        length=127, nullable=False, index=True, default='generic')
    __mapper_args__ = {
        'polymorphic_identity': 'ent',
    }

    # Override default key generator from KeyedMixin
    @classmethod
    def get_key_generator(cls, prefix=None):
        prefix = prefix or cls.__name__.lower()

        def key_generator(ctx):
            ent_type = ctx.current_parameters['ent_type']
            return ':'.join([prefix, ent_type, utils.generate_uuid()])

        return key_generator

    # scoped Ent-Ent relationships.
    children = _orm.relationship(
        'Ent',
        secondary=utils.parent_child,
        primaryjoin=(
            (utils.Node.node_key == utils.parent_child.c.parent_node_key)
            & (utils.Node.node_type == 'ent')
        ),
        secondaryjoin=(
            (utils.Node.node_key == utils.parent_child.c.child_node_key)
            & (utils.Node.node_type == 'ent')
        ),
        backref='parents'
    )

    ancestors = _orm.relationship(
        'Ent',
        secondary=utils.ancestor_descendant,
        primaryjoin=(
            (utils.Node.node_key == (utils.ancestor_descendant.c
                                     .ancestor_node_key))
            & (utils.Node.node_type == 'ent')
        ),
        secondaryjoin=(
            (utils.Node.node_key == (utils.ancestor_descendant.c
                                     .descendant_node_key))
            & (utils.Node.node_type == 'ent')
        ),
        backref='descendants'
    )
