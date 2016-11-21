import psutil

class NoSuchProcess(Exception): pass

def get_process_state(pid=None):
    try:
        process_state = psutil.Process(pid=pid).as_dict()
        process_state['running'] = process_state['status'] not in [
            psutil.STATUS_STOPPED, psutil.STATUS_DEAD, psutil.STATUS_ZOMBIE]
        return process_state
    except psutil.NoSuchProcess as e:
        raise NoSuchProcess(e)
