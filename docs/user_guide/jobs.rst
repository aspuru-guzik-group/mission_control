Jobs
====

==============
What is a job?
==============
In MissionControl, a job represents an standalone computation. For example, a
job that you want to run on a cluster computing resource.

============
Job Lifecyle
============
In general the lifecycle of a job is like this:

#. Create an abstract recipe for job. For example, a dictionary that
   specifies a job type and parameters:
   ::

     {
       'job_type': 'taco',
       'job_params': {
         'wrapper': 'soft_corn',
         'fillings': ['carne', 'salsa', 'lechuga']
       }
     }

#. Generate a runnable artifact that executes computations based on the
   recipe. For example, a directory that contains inputs and a set of commands:
   ::
     taco_job_dir/
       inputs/
         fillings.txt
         wrapper.txt
       entrypoint.sh

     <contents of entrypoint.sh>
     #!/bin/bash
     /bin/taco_maker --fillings="inputs/fillings.txt" wrapper="inputs/wrapper.txt"

#. Run the artifact, usually on a computing cluster. For example:
   ::

     slurm sbatch my_taco_job_dir/entrypoint.sh

#. Store the executed artifact. For example:
   ::

     mv my_taco_job_dir completed_jobs_archive/

#. Parse the executed artifact and store the results in a db. For example:
   ::

     /bin/taco_parser completed_jobs_archive/my_taco_job_dir | /bin/taco_client upload

MissionControl defines conventions and commands to help people manage the
job lifecycle.

job_type
  A string representing a type of job. The job_type is typically used to
  dispatch to specific handlers for a given job. For example, you could specify
  a job_type of 'my_modules.my_taco_module' to dispatch to a python module that
  builds directories for taco jobs.

job_parameters
  A dictionary of parameters. For example, for our taco job, we could specify
  a set of taco parameters, like this:
  ::

    {
      'wrapper': 'soft_corn',
      'fillings': ['carne', 'salsa', 'lechuga']
    }

.. toctree::
   :glob:

   jobs/*
