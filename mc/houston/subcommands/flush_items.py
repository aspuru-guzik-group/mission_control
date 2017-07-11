def flush_mc_db(self, args=None, kwargs=None, unparsed_args=None):
    mc_dao = self._get_mc_dao()
    mc_dao.flush_mc_db()

def flush_flows(self, args=None, kwargs=None, unparsed_args=None):
    mc_dao = self._get_mc_dao()
    mc_dao.delete_items(item_type='Flow')

def flush_jobs(self, args=None, kwargs=None, unparsed_args=None):
    mc_dao = self._get_mc_dao()
    mc_dao.delete_items(item_type='Job')
