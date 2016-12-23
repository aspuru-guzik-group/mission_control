from .base import BaseNode

class A2G2_DAO_Node(BaseNode):
    def tick(self, *args, ctx=None, **kwargs):
        try:
            a2g2_dao = ctx['a2g2_dao']
            query = self.data['input']['query']
            self.data['output'] = {'mol': a2g2_dao.query(query)}
            self.status = 'COMPLETED'
        except Exception as e:
            self.status = 'FAILED'
            self.logger.exception(e)
            self.data['error'] = e
