Jobs
====

=================
Concept Checklist
=================
After you read this section you should be familiar with these concepts:

#. [ ] job lifecycle
#. [ ] job dicts
#. [ ] building job dirs
#. [ ] job modules
#. [ ] running jobs in different environments
#. [ ] parsing job dirs
#. [ ] best practices for working with jobs

==============
What Is a Job?
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

The MissionControl :mod:`mc.utils.job_modules.job_dir_builder`
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
small utility module that is included in MissionControl:
:mod:`mc.utils.testing.echo_job_module` .

It contains a function :mod:`mc.utils.testing.echo_job_module.build_work_dir`
which defines how to build a work_dir.

By convention, MissionControl will look for a function named 'build_work_dir'
in python module that has the same name as the job_type. This function
will receive the job_params and an output_dir as kwargs.

You can also specify a specific builder when you call the 'build_job_dir'
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


============
Running Jobs
============

HERE!!! ADD NOTE ABOUT RUNNING EXTERNALLY, PROVIDING PREBAKED

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Running Jobs in Different Environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Our echo job from the example above is simple and should run the same in any
environment.

But what if want to run jobs that do need special configurations, depending on
the environment in which they run?

For example, what if we want to run job that requires a specific version of a
quantum chemistry library? What if we want to run this job on two different
clusters, cluster X and cluster Y?

There are a few strategies we can use to define environment-specific
configurations.

Config Strategy A: Builder Per Environment
------------------------------------------
In this strategy, we write a builder for each environment in which we expect to
run our job.

For example, our code might look something like this:
  ::

    # <chem_builder_a.py>
    def build_work_dir_for_cluster_x(...):
      # define configs for cluster X
      chem_lib_executable = '/cluster/x/software/my_chem_lib-1.0.1'
      entrypoint = _write_entrypoint(chem_lib_executable)
      return {'entrypoint': entrypoint}
      ...

    def build_work_dir_cluster_x(...): ...
      chem_lib_executable = '/cluster/y/bin/my_chem_lib-1.0.1'
      entrypoint = _write_entrypoint(chem_lib_executable)
      return {'entrypoint': entrypoint}
      ...

    def write_entrypoint(chem_lib_executable):
        entrypoint_content = textwrap.dedent(
          '''
          #!/bin/bash
          CHEM_LIB_EXE="{chem_lib_executable}"
          $CHEM_LIB_EXE my_chem_command
          '''
        ).format(chem_lib_executable=chem_lib_excutable)
        ...

And then when we build our job directories, we just specify which builder to
use:

  ::
    import my_chem_builder_a
    # for cluster x
    houston.run_command(
      'build_job_dir',
      job_dict={...},
      build_work_dir_fn=my_chem_builder_a.build_work_dir_for_cluster_x
    )

    # for cluster y
    houston.run_command(
      'build_job_dir',
      job_dict={...},
      build_work_dir_fn=my_chem_builder_a.build_work_dir_for_cluster_y
    )

Pros
~~~~
#. It's often easier for new users of our code to add new code. "I just copy
   from the previous example!"


Cons
~~~~
#. Maximizing Cluster Use:
   #. We have to know where our job will run at the time we build it. This
      means we would have to check cluster availability at job build time,
      rather than at job run time.
   #. We can't make batches of heterogenous jobs ahead of run time, because
      we would have to check that all the jobs have been built for the same
      cluster.
#. Maintenance:
   #. If our '_generate_common_content' function signature changes,
      we will have to find all the places where is called.
   #. If another type of job uses the same chemstry library, we will have to
      repeat our configurations in the builder for that type of job.
#. Testing: we have to test each of our builders.


Config Strategy B: One Builder + Config Spec
---------------------------------------------
Another strategy is to define one builder, and output a 'config spec' along
with the job_dir. The config spec describes what things this job needs to run.

For example:
  ::

    # <my_chem_builder_b.py>
    def build_work_dir(...):
        # define config spec
        config_spec = {
            'chem_lib_executable': {
                'required': True,
                'env_var': 'CHEM_LIB_EXE'
            }
        }
        return {'entrypoint': entrypoint, 'config_spec': config_spec}

    def write_entrypoint():
        entrypoint_content = textwrap.dedent(
          '''
          #!/bin/bash
          $CHEM_LIB_EXE my_chem_command
          '''
        )
        ...

Pros
~~~~
#. Maintenance: all our logic is one place, so it's easier to maintain.
#. We don't have to know where our job will be run when we build it. So we could
   send it to any cluster that has available resources. And we can batch
   together any collection of jobs.
#. Testing: we only have one builder to test.

Cons
~~~~
#. Whatever runs our job now bears the responsibility for fulfilling the config
   spec requirements.
#. It can be harder for novice users to understand how configs get set.


A vs. B: Which One to Choose?
-----------------------------
In general, the MissionControl authors recommend strategy B. The advantages in
testing and cluster use make up for the slightly higher barrier-to-entry for
job module writers.


================
Parsing Job Dirs
================
Often we want to extract data from executed job dirs.

The MissionControl :mod:`mc.utils.job_modules.job_dir_parser` module
defines a python class that helps us parses a given job directory.

Typically we don't interact directly with the JobDirParser. Instead we use the
`houston` command runner to call the JobDirParser. For example:

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


===========================================
Recommended Practices for Working with Jobs
===========================================
#. Write small functions in your modules.

   This will make your job modules easier to test and understand.

#. Use constants.py files in your modules.

   If your parsers and builders need to refer to common paths or settings, put
   the settings in a constants.py module that both your parsers and builders
   can access. Then, if you need to change these settings, you only need to
   change them in one place.

#. Write tests for your job modules.

#. Define a runner with prebaked outputs.

   This will make your job modules easier to test, both individually and in
   flows.

#. Use the 'One Builder + Config Spec' strategy to specify requirements that
   vary across environments.
