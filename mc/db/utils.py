import json
import uuid
import time

import sqlalchemy as _sqla
import sqlalchemy.orm as _orm
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm.collections import attribute_mapped_collection

from .base import class_registry, Base
from . import constants


def _time(): return time.time()


def _created_or_time(context):
    return context.current_parameters.get('created') or _time()


def generate_int_column(*args, length=None, **kwargs):
    return _sqla.Column(*args, _sqla.types.Integer, **kwargs)


def generate_str_column(*args, length=None, **kwargs):
    return _sqla.Column(*args, _sqla.types.String(length=length), **kwargs)


def generate_boolean_column(*args, **kwargs):
    return _sqla.Column(*args, _sqla.types.Boolean(),
                        **{'default': False, **kwargs})


def generate_json_column(*args, **kwargs):
    return _sqla.Column(*args, JSONEncodedDict(),
                        **{'default': {}, 'nullable': True, **kwargs})


class JSONEncodedDict(_sqla.TypeDecorator):
    impl = _sqla.VARCHAR

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = json.dumps(value, sort_keys=True)
        return value

    def process_result_value(self, value, dialect):
        if value is not None:
            value = json.loads(value)
        return value


def generate_status_column(*args,  **kwargs):
    return _sqla.Column(*args, _sqla.types.String(length=31),
                        **{'default': 'PENDING', 'nullable': True, **kwargs})


def generate_uuid(): return str(uuid.uuid4())


class KeyedMixin(object):
    @declared_attr
    def key(cls):
        return generate_str_column(
            'key',
            length=constants.KEY_LENGTH,
            **{
                'primary_key': True,
                'default': cls.get_key_generator()
            }
        )

    @classmethod
    def get_key_generator(cls, prefix=None):
        prefix = prefix or cls.__name__.lower()
        def key_generator(ctx):  # noqa
            return ':'.join([prefix, generate_uuid()])
        return key_generator

    def __repr__(self):
        return '<{cls}|ID:{id}|KEY:{key}>'.format(
            cls=self.__class__.__name__, id=id(self), key=self.key)


parent_child = _sqla.Table(
    "parent_child", Base.metadata,
    _sqla.Column("parent_node_key", _sqla.String,
                 _sqla.ForeignKey("node.node_key"), primary_key=True),
    _sqla.Column("child_node_key", _sqla.String,
                 _sqla.ForeignKey("node.node_key"), primary_key=True)
)


ancestor_descendant = _sqla.Table(
    "ancestor_descendant", Base.metadata,
    _sqla.Column("ancestor_node_key", _sqla.String,
                 _sqla.ForeignKey("node.node_key"), primary_key=True),
    _sqla.Column("descendant_node_key", _sqla.String,
                 _sqla.ForeignKey("node.node_key"), primary_key=True)
)


class Node(Base):
    node_key = generate_str_column(
        length=constants.KEY_LENGTH,
        **{
            'primary_key': True,
            'default': KeyedMixin.get_key_generator(prefix='node')
        }
    )
    node_type = generate_str_column('node_type')
    child_nodes = _orm.relationship(
        'Node',
        secondary=parent_child,
        primaryjoin=(node_key == parent_child.c.parent_node_key),
        secondaryjoin=(node_key == parent_child.c.child_node_key),
        backref='parent_nodes'
    )
    descendant_nodes = _orm.relationship(
        'Node',
        secondary=ancestor_descendant,
        primaryjoin=(node_key == ancestor_descendant.c.ancestor_node_key),
        secondaryjoin=(node_key == ancestor_descendant.c.descendant_node_key),
        backref='ancestor_nodes'
    )

    __mapper_args__ = {
        'polymorphic_identity': 'node',
        'polymorphic_on': 'node_type',
    }

    def child_nodes_of_type(self, type_=None):
        return [child_node for child_node in self.child_nodes
                if child_node.__class__.__name__ == type_]

    def parent_nodes_of_type(self, type_=None):
        return [parent_node for parent_node in self.parent_nodes
                if parent_node.__class__.__name__ == type_]


class NodeSubclassMixin(object):
    @declared_attr
    def node_ref(cls):
        return _sqla.Column(
            _sqla.String(length=constants.KEY_LENGTH),
            _sqla.ForeignKey('node.node_key'),
        )

    # These add_*_by_key  methods are here for update actions, when we
    # only have a serialized key for a relationship instead of an instance.
    def add_parents_by_key(self, session=None, keys=None):
        with session.begin_nested():
            for key in keys:
                parent = self._get_instance_by_key(session=session, key=key)
                self.parents.append(parent)

    def _get_instance_by_key(self, session=None, key=None):
        Cls = class_registry[key.split(':')[0].title()]
        return session.query(Cls).filter_by(key=key).first()

    def add_ancestors_by_key(self, session=None, keys=None):
        with session.begin_nested():
            for key in keys:
                ancestor = self._get_instance_by_key(session=session, key=key)
                self.ancestors.append(ancestor)


class TimestampMixin(object):
    created = _sqla.Column('created', _sqla.types.Float(),
                           default=_time)
    modified = _sqla.Column('modified', _sqla.types.Float(),
                            default=_created_or_time,
                            onupdate=_time)


class PolymorphicVerticalProperty(object):
    """A key/value pair with polymorphic value storage. """

    def __init__(self, key=None, value=None):
        self.key = key
        self.value = value

    @hybrid_property
    def value(self):
        fieldname, discriminator = self.type_map[self.type]
        if fieldname is None:
            return None
        else:
            return getattr(self, fieldname)

    @value.setter
    def value(self, value):
        py_type = type(value)
        if py_type not in self.type_map:
            fieldname, discriminator = self.type_map['DEFAULT']
        else:
            fieldname, discriminator = self.type_map[py_type]
        self.type = discriminator
        if fieldname is not None:
            setattr(self, fieldname, value)

    @value.deleter
    def value(self): self._set_value(None)

    @value.expression
    def value(self):
        pairs = set(self.type_map.values())
        whens = [
            (
                _sqla.literal_column("'%s'" % discriminator),
                getattr(self, attribute)
            ) for attribute, discriminator in pairs
            if attribute is not None
        ]
        return _sqla.case(whens, self.type, _sqla.null())


@_sqla.event.listens_for(PolymorphicVerticalProperty, "mapper_configured",
                         propagate=True)
def on_new_class(mapper, cls_):
    """Look for Column objects with type info in them, and work up
    a lookup table."""
    info_dict = {}
    info_dict[type(None)] = (None, 'none')
    info_dict['none'] = (None, 'none')

    for k in mapper.c.keys():
        col = mapper.c[k]
        if 'type' in col.info:
            python_type, discriminator = col.info['type']
            info_dict[python_type] = (k, discriminator)
            info_dict[discriminator] = (k, discriminator)
    cls_.type_map = info_dict


class Prop(PolymorphicVerticalProperty):
    key = _sqla.Column(_sqla.String(length=255), primary_key=True)
    type = _sqla.Column(_sqla.String(length=16))
    int_value = _sqla.Column(_sqla.Integer, info={'type': (int, 'integer')})
    float_value = _sqla.Column(_sqla.Float, info={'type': (float, 'float')})
    char_value = _sqla.Column(_sqla.UnicodeText,
                              info={'type': (str, 'string')})
    bool_value = _sqla.Column(_sqla.Boolean, info={'type': (bool, 'boolean')})
    json_value = generate_json_column(
        'json', info={'type': ('DEFAULT', 'json')})

    def __repr__(self):
        return '<{cls} {parent_key}...|{key}={value!r}>'.format(
            cls=self.__class__.__name__,
            parent_key=self.parent_key[:9],
            key=self.key,
            value=self.value
        )


class PropsMixin(object):
    """Creates a new Prop class for each parent."""
    @declared_attr
    def props_set(cls):
        cls.Prop = type(
            "%sProp" % cls.__name__,
            (Prop, Base,),
            {
                '__tablename__': "%s_prop" % cls.__tablename__,
                'parent_key': _sqla.Column(
                    _sqla.String(length=constants.KEY_LENGTH),
                    _sqla.ForeignKey("%s.key" % cls.__tablename__),
                    primary_key=True
                ),
                'parent': _orm.relationship(cls)
            }
        )
        return _orm.relationship(
            cls.Prop,
            collection_class=attribute_mapped_collection('key'),
            cascade='all, delete-orphan'
        )

    @declared_attr
    def props(cls):
        return association_proxy(
            'props_set', 'value',
            creator=(lambda key, value: cls.Prop(key=key, value=value))
        )


class Tag(KeyedMixin, TimestampMixin):
    name = _sqla.Column(_sqla.String(length=255), primary_key=True)


class TagsMixin(object):
    """Creates a new Tag class for each parent."""
    @declared_attr
    def tags_set(cls):
        cls.Tag = type(
            "%sTag" % cls.__name__,
            (Tag, Base,),
            {
                '__tablename__': "%s_tag" % cls.__tablename__,
                'parent_key': _sqla.Column(
                    _sqla.String(length=constants.KEY_LENGTH),
                    _sqla.ForeignKey("%s.key" % cls.__tablename__)
                ),
                'parent': _orm.relationship(cls)
            }
        )
        return _orm.relationship(cls.Tag, collection_class=set)

    @declared_attr
    def tags(cls):
        return association_proxy(
            'tags_set', 'name',
            creator=(lambda name: cls.Tag(name=name))
        )


class ToDictMixin(object):
    def to_dict(self, keys_to_exclude=None):
        keys_to_exclude = keys_to_exclude or set()
        dict_ = {}
        for c in self.__table__.columns:
            if c.name.startswith("_") or c in keys_to_exclude:
                continue
            dict_[c.name] = getattr(self, c.name)
        return dict_


common_supers = [KeyedMixin, TimestampMixin, PropsMixin, TagsMixin,
                 ToDictMixin, NodeSubclassMixin, Node]
