import os

class DjDao(object):
    def __init__(self, DATABASES=None, models=None, serializers=None,
                 queue_utils=None):
        DATABASES = DATABASES or {
            'default': {
                'ENGINE': 'django.db.backends.sqlite3',
                'NAME': ':memory:'
            }
        }
        import django
        from django.conf import settings
        settings.configure(DATABASES=DATABASES, INSTALLED_APPS=['mc.missions'])
        django.setup()
        from django.core.management import call_command
        with open(os.devnull, 'w') as f: call_command('migrate', stdout=f)

        if not models: from mc.missions import models
        self.models = models
        models.Flow(label='foo').save()
        print(models.Flow.objects.all())
        
        if not serializers: from mc.missions import serializers
        self.serializers = serializers

        if not queue_utils: from mc.missions.utils import queue_utils
        self.queue_utils = queue_utils

if __name__ == '__main__':
    dao = DjDao()
