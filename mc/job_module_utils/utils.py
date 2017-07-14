import json
import os

from . import constants

def load_from_jobdir_file(jobdir=None, file_name=None):
    path = os.path.join(jobdir, file_name)
    with open(path) as f: return json.load(f)

def load_job_spec_from_jobdir(jobdir=None):
    return load_from_jobdir_file(
        jobdir=jobdir, file_name=constants.JOB_META_FILE_NAMES['job_spec'])

def load_job_from_jobdir(jobdir=None, file_name=None):
    return load_from_jobdir_file(
        jobdir=jobdir, file_name=constants.JOB_META_FILE_NAMES['job'])
