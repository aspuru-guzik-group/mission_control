import os


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

MC_DB_URI = 'sqlite:///' + os.path.join(ROOT_DIR, 'mc.db.sqlite')
