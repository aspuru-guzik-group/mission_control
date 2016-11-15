import random
from fireworks import Firework, FWorker, LaunchPad, ScriptTask
from fireworks.core.rocket_launcher import launch_rocket

from mc.celery import app as celery_app

LPAD_CFG_FILE = '/mc/fireworks/launchpad_config.yaml'
TAGS = ['a', 'b', 'c']

@celery_app.task
def add(x, y):
    return x + y

@celery_app.task
def add_firework():
    lpad = LaunchPad.from_file(LPAD_CFG_FILE)
    tag = random.choice(TAGS)
    script_str = ('date | tee -a {tag_file}; echo {tag}'
                  ' | tee -a {tag_file}').format(
                      tag=tag, tag_file='/tmp/%s' % tag)
    script_task = ScriptTask.from_str(script_str)
    fw = Firework([script_task])
    fw.spec['_category'] = tag
    lpad.add_wf(fw)

@celery_app.task
def run_fworker():
    lpad = LaunchPad.from_file(LPAD_CFG_FILE)
    tag = random.choice(TAGS)
    fworker = FWorker(category=tag)
    launch_rocket(lpad, fworker)

celery_app.conf.beat_schedule['add_firework_every_5s'] =  {
    'task': add_firework.name,
    'schedule': 5.0,
}

celery_app.conf.beat_schedule['run_rlaunch_every_6s'] =  {
    'task': run_fworker.name,
    'schedule': 6.0,
}
