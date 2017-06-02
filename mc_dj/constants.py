class Statuses(object):
    PENDING = 'PENDING'
    RUNNING = 'RUNNING'
    COMPLETED = 'COMPLETED'
    FAILED = 'FAILED'

    tickable_statuses = [PENDING, RUNNING]
