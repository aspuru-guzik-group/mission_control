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

#. Define a job in temrs of an abstract recipe. For example, a dictionary that
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

=============
Defining Jobs
=============
MissionControl defines jobs as dictionaries with a few primary components:

job_type
  A string representing a type of job, typically used to
  dispatch to specific handlers during the job lifecyle.
  
  For example, you could specify a job with a job_type of
  'my_modules.my_taco_module' to tell MissionControl to use the code in a module
  named 'my_modules.my_taco_module' when it tries to build and parse the job.

job_parameters
  A dict of parameters. For example, for our taco job, we could specify
  a set of taco parameters, like this:
  ::

    {
      'wrapper': 'soft_corn',
      'fillings': ['carne', 'salsa', 'lechuga']
    }

key
  Often you will also want a key to uniquely identify a job. For example, a
  UUID. This key is useful for tracking job statuses, and for avoiding name
  collisions when you are creating and moving multiple job directories.


=============
Building Jobs
=============
Now we have a convention for defining jobs in terms of job_dicts. But how do we
build a runnable artifact for a job_dict?

The MissionControl :mod:`mc.utils.job_modules.job_dir_builder.JobDirBuilder`
module defines a python class that builds a job directory for a given job_dict.

In general JobDirBuilder builds a directory that contains:

#. A work_dir that contains scripts and inputs needed to run the job.
#. An entrypoint script to run the work_dir
#. job metadata files for tracking the job.
#. metadata file(s) that specify what resources or environment variables a job
   needs to run.

Typically we don't interact directly with the JobDirBuilder. Instead we use the
`houston` command runner to call the JobDirBuilder. For example:

.. testcode::

  from mc.houston import Houston
  houston = Houston.minimal()
  import tempfile
  scratch_dir = tempfile.mkdtemp()
  from pathlib import Path
  job_dir_path = Path(scratch_dir, 'my_job_dir')
  build_result = houston.run_command(
     'build_job_dir',
     job_dict={
         'key': 'my_job_key',
         'job_type': 'mc.utils.testing.echo_job_module',
         'job_params': {'message': 'Tacos are delicious.'},
     },
     output_dir=str(job_dir_path)
  )
  built_job_dir = build_result['job_dir']
  print(Path(built_job_dir).name)

The above code creates a Houston instance, and then runs the command
'build_job_dir'. The output of this command is the path to a job dir:

.. testoutput::

  my_job_dir

Let's look at what is inside the job_dir:

.. testcode::

  job_dir_items = [
     str(item_path.relative_to(job_dir_path))
     for item_path in job_dir_path.glob('**/*')
  ]
  print("\n".join(sorted(job_dir_items)))

Expected output:

.. testoutput::

  JOBMAN__JOB_SPEC.json
  MC__JOB_KEY
  MC__JOB_META.json
  entrypoint.sh
  work_dir
  work_dir/entrypoint.sh

We see a list of metadata files and the work_dir .

===========
Job Modules
===========
How did MissionControl know how to build this job dir? The key is the
'job_type' component of the job_dict.

In the example above, we specified a job_type of
'mc.utils.test.echo_job_module'. When we ran the 'build_job_dir' command,
MissionControl looked at the job_type component, and saw that it should try to
dispatch to a module named 'mc.utils.testing.echo_job_module'. This module is a
small utility module that is part of MissionControl
(:mod:`mc.utils.testing.echo_job_module`) .

It contains a function :mod:`mc.utils.testing.echo_job_module.build_work_dir`
which defines how to build a work_dir.

By convention, MissionControl will look for a function named 'build_work_dir'
in python module that has the same name as the job_type. But you can also
specify a specific builder when you call the 'build_job_dir'
command. For example:

.. testcode ::

  def my_build_work_dir(*args, params=None, output_dir=None, **kwargs):
      from pathlib import Path
      Path(output_dir).mkdir(parents=True, exist_ok=True)
      entrypoint_name = 'entrypoint.sh'
      import textwrap
      entrypoint_content = textwrap.dedent(
          '''
          #!/bin/bash
          echo "from my_work_dir_builder"
          echo {message}
          '''
      ).lstrip().format(message=params['message'])
      entrypoint_path = Path(output_dir, entrypoint_name)
      with open(str(entrypoint_path), 'w') as f:
          f.write(entrypoint_content)
      entrypoint_path.chmod(0x775)
      return {'dir': output_dir, 'entrypoint_name': entrypoint_name}

  build_result = houston.run_command(
     'build_job_dir',
     job_dict={
         'key': 'my_job_key',
         'job_type': 'mc.utils.testing.echo_job_module',
         'job_params': {'message': 'Tacos are delicious.'},
     },
     build_work_dir_fn=my_build_work_dir
  )
