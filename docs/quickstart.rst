Quickstart
==========

=======
Install
=======

#. Install via pip:
   ::
     pip install git+https://github.com/aspuru-guzik-group/mission_control.git

#. Verify the install:

   python -m mc.houston.cli sanity_check


===============
Build a Job Dir
===============

.. testcode::

  from mc.houston import Houston
  houston = Houston.minimal()
  build_result = houston.run_command(
     'build_job_dir',
     job_dict={
         'key': 'some_job_key',
         'job_type': 'mc.utils.testing.echo_job_module',
         'job_params': {'message': 'Tacos are delicious.'},
     }
  )
  from pathlib import Path
  job_dir_path = Path(build_result['job_dir'])
  job_dir_items = [
     str(item_path.relative_to(job_dir_path))
     for item_path in job_dir_path.glob('**/*')
  ]
  print("\n".join(sorted(job_dir_items)))
  import subprocess
  entrypoint_command = job_dir_path / 'work_dir' / 'entrypoint.sh'
  job_output = subprocess.check_output([str(entrypoint_command)]).decode()
  print(job_output.strip())

.. testoutput::

  JOBMAN__JOB_SPEC.json
  MC__JOB_KEY
  MC__JOB_META.json
  entrypoint.sh
  work_dir
  work_dir/entrypoint.sh
  Tacos are delicious.

===========
Run a Flow
===========

.. testcode::

  from mc.houston import Houston
  houston = Houston.minimal()
  flow_spec = {
     'tasks': [
         {
             'key': 'task_0',
             'task_type': 'print',
             'task_params': {'msg': 'Hello from task_0'}
         },
         {
             'key': 'task_1',
             'task_type': 'print',
             'task_params': {'msg': 'Hello from task_1'},
             'data': {
                 'my_msg': 'Message for task_2, from task_1.data'
             }
         },
         {
             'key': 'task_2',
             'task_type': 'print',
             'task_params': {'msg': '$ctx.flow.tasks.task_1.data.my_msg'}
         }
     ]
  }
  completed_flow_dict = houston.run_command(
     'tick', tickee='flow',
     flow_spec=flow_spec,
     until_finished=True,
     interval=.01,
  )
  print(completed_flow_dict['status'])

.. testoutput::

  Hello from task_0
  Hello from task_1
  Message for task_2, from task_1.data
  COMPLETED
