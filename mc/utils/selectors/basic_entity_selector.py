import sqlalchemy as _sqla
import sqlalchemy.orm as _orm

from .base_selector import BaseSelector


class BasicEntitySelector(BaseSelector):
    PAGE_SIZE = int(1e3)

    def __init__(self, db=None, **kwargs):
        self.db = db

    @property
    def Ent(self): return self.db.models.Ent

    @property
    def session(self): return self.db.session

    def get_items(self, ent_type=None, having_tags=None, sans_tags=None,
                  sans_requests_w_tags=None, **kwargs):
        q = self.session.query(self.Ent).filter_by(ent_type=ent_type)
        if having_tags:
            q = self._alter_q_per_having_tags(q=q, having_tags=having_tags)
        if sans_tags:
            q = self._alter_q_per_sans_tags(q=q, sans_tags=sans_tags)
        if sans_requests_w_tags:
            q = self._alter_q_per_sans_requests_w_tags(
                q=q, sans_requests_w_tags=sans_requests_w_tags)
        q = q.order_by(self.Ent.modified)
        for entity in q.yield_per(self.PAGE_SIZE):
            yield {'key': entity.key, 'value': entity}

    def _alter_q_per_having_tags(self, q=None, having_tags=None):
        _Tag = self.Ent.Tag
        tag_subq = (
            self.db.session.query(_Tag)
            .filter(_Tag.name.in_(having_tags))
            .group_by(_Tag.parent_key)
            .having(_sqla.func.count(_Tag.key) == len(having_tags))
        )
        return q.join(tag_subq.subquery())

    def _alter_q_per_sans_tags(self, q=None, sans_tags=None):
        return q.filter(
            ~(self.Ent.tags_set.any(self.Ent.Tag.name.in_(sans_tags)))
        )

    def _alter_q_per_sans_requests_w_tags(
        self, q=None, sans_requests_w_tags=None):  # noqa
        _Request = _orm.aliased(self.db.models.Request)
        _instance_key = _Request.instance_key.label('instance_key')
        requests_subq = (
            self.db.session.query(_instance_key)
            .filter(_Request.request_tag.in_(sans_requests_w_tags))
            .subquery()
        )
        q = (
            q.outerjoin(
                requests_subq,
                (self.Ent.key == requests_subq.c[_instance_key.name])
            )
            .group_by(self.Ent)
            .having(_sqla.func.count(requests_subq.c[_instance_key.name]) == 0)
        )
        return q
