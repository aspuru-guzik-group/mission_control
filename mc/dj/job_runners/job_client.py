from enum import Enum

class JobClient(object):
    class Statuses(Enum):
        Completed = 1
        Transferred = 2

    def fetch_jobs(self):
        pass

    def claim_jobs(self, uuids=None):
        pass

    def update_jobs(self, updates_by_uuid=None):
        pass

