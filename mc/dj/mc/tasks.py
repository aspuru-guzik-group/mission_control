from mc.celery import app as celery_app

@celery_app.task
def add(x, y):
    return x + y

@celery_app.task
def add_firework():
    from fireworks import Firework, LaunchPad, ScriptTask
    lpad_cfg_file = '/mc/fireworks/launchpad_config.yaml'
    lpad = LaunchPad.from_file(lpad_cfg_file)
    script_task = ScriptTask.from_str('date >> /tmp/beef.txt;'
                                      'echo "beef" >> /tmp/beef.txt"')
    fw = Firework([script_task])
    lpad.add_wf(fw)
