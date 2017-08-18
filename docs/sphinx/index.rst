.. MissionControl documentation master file, created by
   sphinx-quickstart on Tue Jun 27 15:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: Introduction to MissionControl (workflow software)

==============
MissionControl
==============
MissionControl: A library for jobs and workflows.

========
Features
========

#. Generate jobs for computing clusters: write job modules that
   build and parse submissions for computing clusters.
#. WorkFlows: create dynamic workflows that run tasks.
#. Flexible Database Backend: use any SqlAlchemy-compatible database to track
   jobs and workflows. Use server-based databases like Postgresql or MySql,
   or store your data in a file or in memory with Sqlite.
#. Store and Query Domain Data: Use an EAV data model to store and query a wide
   range of records.
#. Modular: Use any combination of MissionControl components you want.

==========
Quickstart
==========

-------
Install
-------

#. Install via pip:
   ::
     pip install git+https://github.com/aspuru-guzik-group/mission_control.git

#. Verify the install:

   python -m mc.houston.cli sanity_check


---------------
Build a Job Dir
---------------

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
  import subprocess
  entrypoint_command = job_dir_path / 'work_dir' / 'entrypoint.sh'
  job_output = subprocess.check_output([str(entrypoint_command)]).decode()
  print("= JOB OUTPUT =")
  print(job_output)
  print("= JOB DIR CONTENTS =")
  job_dir_items = [
     str(item_path.relative_to(job_dir_path))
     for item_path in job_dir_path.glob('**/*')
  ]
  print("\n".join(sorted(job_dir_items)))

Output:

.. testoutput::

  = JOB OUTPUT =
  Tacos are delicious.

  = JOB DIR CONTENTS =
  MC__JOB_KEY
  MC__JOB_META.json
  entrypoint.sh
  job_spec.json
  work_dir
  work_dir/echo_job.out
  work_dir/entrypoint.sh

-----------
Run a Flow
-----------

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

---------------------------
Store and Query Domain Data
---------------------------

.. testcode::

  from mc.houston import Houston
  houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
  houston.db.ensure_tables()

  molecules = []
  Ent = houston.db.models.Ent
  for i in range(1, 4):
      molecule = Ent(
          ent_type='molecule',
          props={'num_atoms': i}
      )
      molecules.append(molecule)
  houston.db.session.add_all(molecules)
  houston.db.session.commit()

  molecules_w_multiple_atoms = (
      houston.db.session.query(Ent)
      .join(Ent.Prop, aliased=True, from_joinpoint=True)
      .filter(Ent.Prop.key == 'num_atoms')
      .filter(Ent.Prop.value > 1)
      .order_by(Ent.Prop.value)
      .all()
  )
  for molecule in molecules_w_multiple_atoms:
     print('num_atoms:', molecule.props['num_atoms'])

.. testoutput ::

  num_atoms: 2
  num_atoms: 3

===============
Getting Started
===============

.. toctree::
   :maxdepth: 1

   user_guide
   examples


=================
API Documentation
=================

.. toctree::
   :maxdepth: 1

   api

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
