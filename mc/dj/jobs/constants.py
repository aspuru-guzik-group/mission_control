import enum

class JobStatuses(enum.Enum):
    PENDING = {'label': 'pending'}
    CLAIMED = {'label': 'claimed'}
    COMPLETED = {'label': 'completed'}
