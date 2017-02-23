import enum


class JobStatuses(enum.Enum):
    PENDING = {'label': 'pending'}
    RUNNING = {'label': 'running'}
    COMPLETED = {'label': 'completed'}
    FAILED = {'label': 'failed'}
