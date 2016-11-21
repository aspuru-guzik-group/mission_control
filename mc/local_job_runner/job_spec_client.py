from enum import Enum

class JobSpecClient(object):
    class Statuses(Enum):
        COMPLETED = 1
        TRANSFERRED = 2

    def fetch_job_specs(self):
        pass

    def claim_job_specs(self):
        pass

    def update_job_specs(self):
        pass

