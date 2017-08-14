Quickstart
==========

=======
Install
=======

#. Install via pip::

   pip install git+https://github.com/aspuru-guzik-group/mission_control.git

#. Verify the install:

   python -m mc.houston.cli sanity_check


===========
Build a Job
===========

.. testcode::

  from mc.houston import Houston
  houston = Houston(ensure_db=False, ensure_job_dirs=False)
  job_dir = houston.run_command(
     'build_job_dir',
     job_dict={
         'key': 'some_job_key',
         'job_type': 'mc.utils.testing.echo_job',
         'job_params': {'message': 'Tacos are delicious.'},
     }
  )
  import glob
  glob.glob(job_dir + '/**/*', recursive=True)


===========
Run a Flow
===========

.. testcode::

  from mc.houston import Houston
  houston = Houston(ensure_db=False, ensure_job_dirs=False)
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
  completed_flow_spec = houston.run_command(
     'tick', tickee='flow',
     flow_spec=flow_spec,
     until_completed=True,
  )
