Jobs
====

A MissionControl job represents a standalone computation. For example, a
single run of computational chemistry model.

============
Job Lifecyle
============
In general the lifecycle of a job is like this:

#. Define a job in terms of an abstract recipe. Example: a dictionary that
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
   recipe. Example: a directory that contains inputs and a set of commands:
   ::

     taco_job_dir/
       inputs/
         fillings.txt
         wrapper.txt
       entrypoint.sh

     <contents of entrypoint.sh>
     #!/bin/bash
     /bin/taco_maker --fillings="inputs/fillings.txt" wrapper="inputs/wrapper.txt"

#. Run the artifact locally or on a computing cluster. Example:
   ::

     slurm sbatch my_taco_job_dir/entrypoint.sh

#. Store the executed artifact. Example:
   ::

     mv my_taco_job_dir completed_jobs_archive/

#. Parse the executed artifact and store the results in a db. Example:
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
  'my_modules.my_taco_module'. Then MissionControl tries to build and parse the
  job, it would know to use code in a module named 'my_modules.my_taco_module'.

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
  collisions when you create and move multiple jobs.


=============
Building Jobs
=============
Now we have a convention for defining jobs in terms of job_dicts. But how do we
build a runnable artifact for a job_dict?

The MissionControl :mod:`mc.utils.job_modules.job_dir_builder`
module defines a python class that builds a job directory for a given job_dict.

JobDirBuilder builds a directory that contains:

#. A work_dir that contains scripts and inputs needed to run the job.
#. An entrypoint script to run the work_dir
#. job metadata files for tracking the job.
#. metadata file(s) that specify what resources or environment variables a job
   needs to run.

We usually don't interact directly with the JobDirBuilder. Instead we use the
:doc:`houston` command runner to call the JobDirBuilder. For example:

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
'build_job_dir'. The output of the above example is the path to a job dir:

.. testoutput::

  my_job_dir

Let's look at the contents of the job_dir:

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
How did MissionControl know how to build the job dir in the example abeove?
The key is the 'job_type' component of the job_dict.

Here is what happend:

#. In our job_dict we specified a job_type of
'mc.utils.test.echo_job_module'.

#. When we ran the 'build_job_dir' command, MissionControl looked at the
   job_type component, and saw that it should try to dispatch to a module named
   'mc.utils.testing.echo_job_module'. This module is a small utility module
   that is included in MissionControl: :mod:`mc.utils.testing.echo_job_module` . 
   It contains a function :mod:`mc.utils.testing.echo_job_module.build_work_dir`
   which defines how to build a work_dir.

#. By convention, MissionControl loaded a function named 'build_work_dir'
   in the python module that has the same name as the job_type. This function
   received the job_params and an output_dir as kwargs.  

#. MissionControl executed the 'build_work_dir' function

#. The 'build_work_dir' function created the job directory.

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

Once you have built a job how can you run it?

There are several ways to run a job. We outline a few here.

-----------------------
Strategy A: Run Locally
-----------------------
In this strategy, you simply execute bash commands to execute a job. Example:
  ::

    bash my_job_dir/entrypoint.sh

~~~~
Pros
~~~~
#. This is the simplest strategy. It requires no coordination with external
   systems.

~~~~
Cons
~~~~
#. Your jobs are limited by the resources of your local environment. This may
   make it difficult to run large numbers of jobs, or to run jobs that require
   lots of CPUs or special software.

---------------------------------------------------
Strategy B: Run on an External Computation Resource 
---------------------------------------------------
In this strategy, you submit jobs directly to an external resource, such as 
a computing cluster, or an EC2 instance. Example:
  ::

    scp my_job_dir me@cluster_x:my_jobs
    ssh me@cluster_x 'cd my_jobs/my_job_dir && sbatch entrypoint.sh'

~~~~
Pros
~~~~
#. You can take advantage of cluster nodes to run large numbers of jobs.
#. You can use special hardware and software provided by a cluster.

~~~~
Cons
~~~~
#. You need to coordinate submissions with an external resource.
#. You need to monitor job progress.
#. It can take time to transfer job files to and from the external resource.

--------------------------------------
Strategy C: Submit to a meta-scheduler
--------------------------------------

In this strategy, you submit jobs directly to a meta-scheduler, such as
`jobman.` The meta-scheduler runs jobs locally or submits them to external
resources, depending on what resources are available. The meta-scheduler can
also create batches of jobs and optimize file transfers.
  ::

    python -m jobman.cli submit_job_dir my_job_dir

~~~~
Pros
~~~~
#. You can use simultaneously use multiple clusters.
#. You only have to coordinate with the meta-scheduler, as opposed to
   coordinating with several different clusters.
#. You can run jobs more efficiently if the meta-scheduler can make batches
   of jobs.

~~~~
Cons
~~~~
#. You need to coordinate with the meta-scheduler.
#. You need to configure the meta-scheduler.

-------------------------
Which Strategy To Choose?
-------------------------
It depends. In general, choose whichever strategy is easiest to get started
with first.


========================
Configuring How Jobs Run
========================

Often we want to specify how jobs should run. Common run parameters include:

#. How much memory should be allocated for a job?
#. How many cores can a job use?
#. Can a job be batched with other jobs?
#. Which specific executables should a job use?

These parameters are different from job parameters. Job parameters
specify what a job should run ('run model X for three iterations'). Whereas
run parameters specify how a job should run ('run the job using 4 cores').

MissionControl has a convention for specifying run parameters. When
MissionControl builds a job dir, it outputs a 'run_params' metadata file.
This metadata file is a generic JSON file. Various job runners can read this
file and then translate the given run parameters into parameters specific to a
given run environment.

For example, a job runner that runs job dirs in a Slurm environment can
translate parameters for memory and cores into Slurm-specific parameters.

By default MissionControl generates a run_params file that is compatible with
the `Jobman` meta-scheduler. 

-----------------------------------
Environment-Specific Configurations
-----------------------------------

Let's think about our echo job from the examples above. It is simple and should
run the same in any environment.

But what if want to run jobs that need special configurations. What if these
configurations depend on the environment the job runs in?

For example, what if we want to run job that requires a specific version of a
quantum chemistry library? What if we want to run this job on two different
clusters, cluster X and cluster Y?

There are a few strategies we can use to define environment-specific
configurations.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Strategy A: Builder Per Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Strategy B: One Builder + Config Spec
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A vs. B: Which One to Choose?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
:doc:`houston` command runner to call the JobDirParser. For example:

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

  # Here we execute a 'fake' run, using prebaked output.
  # This often a useful strategy for testing parsers.
  houston.run_command(
     'run_job_dir',
     job_dir=built_job_dir,
     fake=True
  )

  parse_results = houston.run_command(
     'parse_job_dir',
     job_dir=built_job_dir,
  )
  print(parse_results)

Expected output:

.. testoutput::

  {'output': 'Tacos are delicious.\n'}


--------------
Parser Outputs
--------------
What should a parser return as outputs? It depends.

MissionControl has no requirements on what a parser must return.

Often what you want a parser to return is some set of update specs that you can
pass to a database client, in order to update records in a database.

The MissionControl `EntityDb` can help you with this type of parsing. See
 `examples.entity_db_parse_and_load`.


===================
Testing Job Modules
===================
The MissionControl authors strongly recommend writing tests for your job module
functions, for several reasons:

#. Debugging is easier in the context of small, isolated tests. Debugging is
   much harder when you have many jobs running in different environments.

#. Tests act as a form of documentation for your code.

#. Tests let you change your code with confidence.

For examples of how to test job module functions, see the source of
:mod:`mc.utils.testing.tests.test_echo_job_module`.

--------------------------
Including Prebaked Outputs
--------------------------
Something that is often helpful for testing jobs is including prebaked outputs
with your job module code. Example: including the output file of a long-running
computational chemistry model alongside your job module code.

Prebaked outputs are helpful for several reasons:

#. prebaked outputs can help with testing parsers.
#. prebaked outputs make it easier to test your jobs in the context of flows
   and pipelines.

For an example of one way to include prebaked outputs, see the source of
:mod:`mc.utils.testing.echo_job_module`.


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
   the context of flows.

#. Use the 'One Builder + Config Spec' strategy to specify requirements that
   vary across environments.

#. Write tests for your job modules.
