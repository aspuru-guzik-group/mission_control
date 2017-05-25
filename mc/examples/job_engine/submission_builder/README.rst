A common task is building a job directory which uses a bash script as its entrypoint. Ideally this entrypoint creates (1) 'checkpoint' files for tracking job completion and failure, (2) standard log files for logging stdout and stderr, and (3) a standard set of directory names for input artifacts and output artifacts.

mc.job_engines provides a small helper module to assist with this task. This module is mc.job_engines.submission_builders.bash .

This module defines a class named BashSubmissionBuilder. This class has a method named build_submission. The build_submission method takes a submission_spec, a cfg object, and an output_dir. It writes an entrypoint file to the output_dir which contains:
1. bash commands to setup checkpoint files and log files
2. bash commands to set environment variables, provided in the submission_spec.
3. A body section, which is provided in the submission_spec.


