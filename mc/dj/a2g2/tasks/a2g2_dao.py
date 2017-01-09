from .base import BaseTask

class A2G2_DAO_Task(BaseTask):
    def tick(self, *args, ctx=None, **kwargs):
        try:
            a2g2_dao = ctx['a2g2_dao']
            if 'query' in self.input:
                query = self.input['query']
                self.output = {'result_set': a2g2_dao.query(query)}
                self.status = 'COMPLETED'
        except Exception as e:
            self.status = 'FAILED'
            self.logger.exception(e)
            self.error = e
