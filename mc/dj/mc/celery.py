import os

from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mc.settings')
app = Celery('mc', broker='redis://redis:6379')
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks()

app.conf.beat_schedule = {}
app.conf.timezone = 'UTC'
