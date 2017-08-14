import sqlalchemy as _sqla
import sqlalchemy.orm as _orm

from .base_selector import BaseSelector


class BasicRequestSelector(BaseSelector):
    PAGE_SIZE = int(1e3)

    def __init__(self, db=None, **kwargs):
        self.db = db

    @property
    def Request(self): return self.db.models.Request

    @property
    def session(self): return self.db.session

    def get_items(self, having_request_tag=None,
                  having_request_type=None, having_status=None,
                  sans_downstream_requests_w_tags=None, **kwargs):
        q = self.session.query(self.Request)
        if having_request_tag:
            q = q.filter(self.Request.request_tag == having_request_tag)
        if having_request_type:
            q = q.filter(self.Request.request_type == having_request_type)
        if having_status:
            q = q.filter(self.Request.status == having_status)
        if sans_downstream_requests_w_tags:
            q = self._alter_q_per_sans_downstream_requests_w_tags(
                q=q,
                sans_downstream_requests_w_tags=sans_downstream_requests_w_tags
            )
        q = q.order_by(self.Request.modified)
        for request in q.yield_per(self.PAGE_SIZE):
            yield {'key': request.key, 'value': request}

    def _alter_q_per_sans_downstream_requests_w_tags(
        self, q=None, sans_downstream_requests_w_tags=None):  # noqa
        _Request = _orm.aliased(self.Request)
        _instance_key = _Request.instance_key.label('instance_key')
        requests_subq = (
            self.db.session.query(_instance_key)
            .filter(_Request.request_tag.in_(sans_downstream_requests_w_tags))
            .subquery()
        )
        q = (
            q.outerjoin(
                requests_subq,
                (self.Request.key == requests_subq.c[_instance_key.name])
            )
            .group_by(self.Request)
            .having(_sqla.func.count(requests_subq.c[_instance_key.name]) == 0)
        )
        return q
