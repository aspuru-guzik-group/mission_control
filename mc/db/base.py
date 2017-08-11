from sqlalchemy.ext.declarative import declared_attr, declarative_base


class Base(object):
    """Base class which provides automated table name."""
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()


Base = declarative_base(cls=Base)
