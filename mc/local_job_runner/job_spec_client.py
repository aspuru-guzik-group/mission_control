from enum import Enum

class JobSpecClient(object):
    class Statuses(Enum):
        COMPLETED = 1
        TRANSFERRED = 2

    def fetch_job_specs(self):
        pass

    def claim_job_specs(self, uuids=None):
        pass

    def update_job_specs(self, updates_by_uuid=None):
        pass

