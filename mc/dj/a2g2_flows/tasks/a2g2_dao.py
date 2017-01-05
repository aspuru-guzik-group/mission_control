from .base import BaseTask

class A2G2_DAO_Task(BaseTask):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)

    def tick(self, *args, ctx=None, **kwargs):
        try:
            a2g2_dao = ctx['a2g2_dao']
            if 'query' in self.data['input']:
                query = self.data['input']['query']
                self.data['output'] = {'result_set': a2g2_dao.query(query)}
                self.status = 'COMPLETED'
            elif 'ingest' in self.data['input']:
                ingest_spec = self.data['input']['ingest']
                a2g2_dao.ingest(ingest_spec=ingest_spec)
                self.status = 'COMPLETED'
        except Exception as e:
            self.status = 'FAILED'
            self.logger.exception(e)
            self.data['error'] = e
