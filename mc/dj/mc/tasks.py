import random
from fireworks import Firework, FWorker, LaunchPad, ScriptTask
from fireworks.core.rocket_launcher import launch_rocket

from mc.celery import app as celery_app
from mc.firetasks.process_job_request import ProcessJobRequestFireTask

LPAD_CFG_FILE = '/mc/fireworks/launchpad_config.yaml'
TAGS = ['a', 'b', 'c']

@celery_app.task
def add(x, y):
    return x + y

def _add_firework(fw):
    lpad = LaunchPad.from_file(LPAD_CFG_FILE)
    lpad.add_wf(fw)

@celery_app.task
def add_firework():
    tag = random.choice(TAGS)
    script_str = ('date | tee -a {tag_file}; echo {tag}'
                  ' | tee -a {tag_file}').format(
                      tag=tag, tag_file='/tmp/%s' % tag)
    script_task = ScriptTask.from_str(script_str)
    fw = Firework([script_task])
    fw.spec['_category'] = tag
    _add_firework(fw)

@celery_app.task
def run_fworker(category=None):
    lpad = LaunchPad.from_file(LPAD_CFG_FILE)
    fworker = FWorker(category=category)
    launch_rocket(lpad, fworker)

def request_molgen():
    ftask = ProcessJobRequestFireTask()
    fwork = Firework([ftask])
    fwork.spec['_category'] = 'job_requests'
    _add_firework(fwork)

celery_app.conf.beat_schedule['job_requests'] =  {
    'task': run_fworker.name,
    'schedule': 10.0,
    'kwargs': {
        'category': 'job_requests',
    }
}

celery_app.conf.beat_schedule['odyssey_job_requests'] =  {
    'task': run_fworker.name,
    'schedule': 12.0,
    'kwargs': {
        'category': 'odyssey_job_requests',
    }
}

celery_app.conf.beat_schedule['odyssey_job_polls'] =  {
    'task': run_fworker.name,
    'schedule': 13.0,
    'kwargs': {
        'category': 'odyssey_job_polls',
    }
}

celery_app.conf.beat_schedule['odyssey_transfers'] =  {
    'task': run_fworker.name,
    'schedule': 15.0,
    'kwargs': {
        'category': 'odyssey_transfers',
    }
}

celery_app.conf.beat_schedule['completed_job_processing'] =  {
    'task': run_fworker.name,
    'schedule': 16.0,
    'kwargs': {
        'category': 'completed_job_processing',
    }
}
