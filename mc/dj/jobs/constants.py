import enum

class JobStatuses(enum.Enum):
    Pending = {'label': 'pending'}
    Claimed = {'label': 'claimed'}
    Completed = {'label': 'completed'}
