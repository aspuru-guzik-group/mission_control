Houston
=======

=================
Concept Checklist
=================
This page covers these concepts:

#. [ ] What is Houston?
#. [ ] Why use Houston?
#. [ ] Running commands with Houston
#. [ ] Houston utils
#. [ ] Configuring Houston


================
What is Houston?
================
:mod:`Houston <mc.houston.houston.Houston>` is a high-level interface to MissionControl.

----------
An Analogy
----------
Think of a car. A car is a system of many components.
There's an engine, brakes, gears, wheels, and many other parts.
But when you drive a car, you don't think about each little part;
it would be too overwhelming! Instead,
you just use the steering wheel and the pedals. These are your interface to
the car.

MissionControl is also a system of many components. There are jobs, flows, 
databases, utility classes, and many other parts. When we use MissionControl,
we don't want to have to know about all of the little parts. We want something
that hides the details. Houston is a set of functions and utilities that lets
us use MissionControl without having to know about all the little parts.
Houston is like a steering wheel for MissionControl.

================
Why use Houston?
================
We use Houston for several reasons:

#. It makes it easier for us to write programs that interact with
   MissionControl.
#. It gives us a way to consolidate configuration in one place.
#. It gives us a framework to extend MissionControl with new commands.


===================
Configuring Houston
===================
Houston takes a 'cfg' object as a parameter. This cfg object can be a dict that
maps config keys to values. Or it can be a python object that has config keys
as attributes.

This allows you to use both normal dicts and python modules as config sources.

--------------
Config as Dict
--------------

.. testcode ::

  from mc.houston import Houston
  my_config = {
     'SOME_CFG_KEY': 'SOME_CFG_VALUE'
  }
  houston = Houston(cfg=my_config)
  print(houston.cfg['SOME_CFG_KEY'])

Expected output:

.. testoutput ::

  SOME_CFG_VALUE


--------------
Config as Object
--------------

.. testcode ::

  import types
  my_config = types.SimpleNamespace()
  my_config.SOME_CFG_KEY = 'SOME_CFG_VALUE'
  houston = Houston(cfg=my_config)
  print(houston.cfg['SOME_CFG_KEY'])

Expected output:

.. testoutput ::

  SOME_CFG_VALUE

------------------
Common Config Keys
------------------
These are some of the common config keys used by Houston.

MC_DB_URI
  The database connection url for the MissionControl database. See
  http://docs.sqlalchemy.org/en/latest/core/engines.html#database-urls for
  valid url formats.
  
FLOW_QUEUE
  Flow queue config. Only needed if you use Houston.utils.flow_runner.

  Example: 
  ::

    {
        'key': 'my_flow_queue',
        'queue_kwargs': {
            'queue_spec': {
                'item_type': 'Flow'
            }
        }
    }

JOB_QUEUE
  Job queue config. Only needed if you use Houston.utils.job_runner.

  Example: 
  ::

    {
        'key': 'my_job_queue',
        'queue_kwargs': {
            'queue_spec': {
                'item_type': 'Job'
            }
        }
    }

USE_LOCKS
  Set to True to enable the use of locks when claiming flows.
  Default: True.

ARTIFACT_HANDLER
  An instance of an artifact handler to use for converting dirs to artifacts.

JOB_DIRS_ROOT
  The root path to use for job dirs.

JOBMAN_CFG
  A jobman cfg specification.


=============================
Running Commands with Houston
=============================
Normally you use Houston to call MissionControl commands.

The typical way to run a command is to use Houston's 'run_command' function.

Example:
  ::

    result = houston.run_command(
       'some_command',
       some_kwarg=...,
       some_other_kwarg=...
    )

Note: there is also a 'call_command' function. That function is intended for
use with command-line calls.

=================
Houston Utilities
=================
If you want to interact directly with MissionControl components you can 
usually access them through Houston's 'utils' attribute. The 'utils' attribute
is an instance of :mod:`mc.houston.houston.utils.HoustonUtils`. You can view
the source of that module to see what things are provided by utils.
