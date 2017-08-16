Flows
=====

==============
What Is a Flow?
==============
A MissionControl Flow represents a set of tasks that run in a specific
sequence.

============
Flow Lifecyle
============
The lifecycle of a flow is like this:

#. Define a flow in terms of an abstract spec. Example: a dictionary that
   specifies a list of tasks:
   ::

     flow_spec = {
         'label': 'example_flow',
         'tasks': [
             {
                 'key': 'task_1',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_1.'},
             },
             {
                 'key': 'task_2',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_2.'},
             },
         ]
     }

#. Convert the flow spec into a python Flow object:
   ::

     flow = flow_spec_to_flow(flow_spec)

#. Tick the flow until either (A) it fails, or (B) all of its tasks finish:
   ::

     while flow.status not in {'FAILED', 'COMPLETED'}:
         tick_flow(flow)

MissionControl defines conventions and commands to help manage the
flow lifecycle.

==============
Defining Flows
==============
Abstract flow specs are defined as dictionaries with a few primary
components:

label
  a human-readable label for the flow

key
  A key that uniquely identifies a flow.
  For example, a UUID. This key is useful for tracking flow statuses, and for
  avoiding name collisions when you work with multiple flows.

cfg
  Configuration parameters for how the flow should run. For example, if one
  task fails should the entire flow fail?

data
  Initial data that the flow should have.

tasks
  A list of task_specs. See :doc:`tasks`.

Example:
   ::

     my_flow_spec = {
         'label': 'example_flow',
         'key': 'some-unique-key-12345',
         'cfg': {
           'fail_fast': False
         },
         'tasks': [
             {
                 'key': 'task_1',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_1.'},
             },
             {
                 'key': 'task_2',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_2.'},
             },
         ]
     }

==============
Building Flows
==============
Now we can define flow_specs. But how do turn our flow_spec into a Flow object?

To convert our flow_spec into a flow object we call
:meth:`mc.flows.flow_engine.FlowEngine.flow_spec_to_flow`.

We get back an instance of :class:`mc.flows.flow.Flow`. This class contains
methods for querying flow tasks and manipulating a flow's underlying
attributes.

=============
Running Flows
=============
Now that we have a Flow object, we can run it.  We use
:class:`FlowEngine<mc.flows.flow_engine.FlowEngine>` to run flows. FlowEngine
is a class which contains methods for ticking flows. It also provides wrappers
for the Flow conversion methods.

.. testcode::

     flow_spec = {
         'label': 'example_flow',
         'tasks': [
             {
                 'key': 'task_1',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_1.'},
             },
             {
                 'key': 'task_2',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_2.'},
             },
         ]
     }
     from mc.flows.flow_engine import FlowEngine
     my_flow_engine = FlowEngine()
     flow = my_flow_engine.flow_spec_to_flow(flow_spec)
     my_flow_engine.tick_flow_until_has_no_pending(flow)
     print("flow.status:", flow.status)

Expected output:

.. testoutput::

   I am task_1.
   I am task_2.
   flow.status: COMPLETED

================
Persistent Flows
================
Often you will want to have several flows which persist over time.

For example, you may want to have a flow runner which runs a loop like this:

#. Retrieves a list of pending flows from a database.
#. Ticks those flows until they fail or have no more pending tasks.
#. Saves the update flows back to the database.

The general lifecycle for storing flows is like this:

#. Serialize a flow into a format suitable for storage.
#. Save the serialized flow to a database.
#. Load serialized flows from the database.
#. Deserialize the serialized flows back to normal flows.

-----------------------------------
Serializing And Deserializing Flows
-----------------------------------

In order to save and load flows, we need to transform flow objects in data
structures which can be stored in a database. The :class:`mc.flows.Flow`
class has class methods for this transformation:

.. testcode:: serialization_test_group

     flow_spec = {
         'label': 'example_flow',
         'tasks': [
             {
                 'key': 'task_1',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_1.'},
             },
             {
                 'key': 'task_2',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_2.'},
             },
         ]
     }
     from mc.flows.flow import Flow
     flow = Flow.from_flow_spec(flow_spec)
     flow_dict = flow.to_flow_dict()
     import json
     jsonified_flow = json.dumps(flow_dict, indent=2, sort_keys=True)
     print("jsonified flow:\n", jsonified_flow)

Expected output:

.. testoutput:: serialization_test_group

    jsonified flow:
     {
      "cfg": {
        "fail_fast": true
      },
      "data": {},
      "depth": 0,
      "graph": {
        "edges": [
          {
            "dest_key": "task_1",
            "src_key": "ROOT"
          },
          {
            "dest_key": "task_2",
            "src_key": "task_1"
          }
        ],
        "tasks": {
          "ROOT": {
            "key": "ROOT",
            "status": "COMPLETED"
          },
          "task_1": {
            "key": "task_1",
            "precursors": [
              "ROOT"
            ],
            "status": "PENDING",
            "task_params": {
              "msg": "I am task_1."
            },
            "task_type": "print"
          },
          "task_2": {
            "key": "task_2",
            "precursors": [
              "task_1"
            ],
            "status": "PENDING",
            "task_params": {
              "msg": "I am task_2."
            },
            "task_type": "print"
          }
        }
      },
      "label": "example_flow",
      "num_tickable_tasks": 1,
      "parent_key": null,
      "status": "PENDING"
    }

Notice how the serialized flow represents the flow's underlying graph.

To deserialize the serialized flow, we can do something like this:

.. testcode:: serialization_test_group

     flow_dict = json.loads(jsonified_flow)
     flow = Flow.from_flow_dict(flow_dict)
     print(type(flow))

.. testoutput:: serialization_test_group

  <class 'mc.flows.flow.Flow'>


------------------------
Saving and Loading Flows
------------------------
MissionControl provides utilities for saving flow_dicts to a database, and
for querying flows.

These utilities are provided by MissionControl's
:doc:`Houston <houston>` utility.

~~~~~~~~~~~~
Saving Flows
~~~~~~~~~~~~

We can save flows using SqlAlchemy actions using Houston's db utility.

.. testcode:: flow_db_test_group

   # Setup houston w/ an in-memory sqlite db.
   from mc.houston import Houston
   my_houston = Houston(
       cfg={
           'MC_DB_URI': 'sqlite://'
       }
   )
   my_houston.db.ensure_tables()

   # Create a flow.
   flow_spec = {
       'label': 'example_flow',
       'tasks': [
           {'task_type': 'print', 'task_params': {'msg': 'I am task_%s.' % i}}
           for i in range(3)
        ]
   }
   from mc.flows.flow import Flow
   flow = Flow.from_flow_spec(flow_spec)

   # Save the flow to the db.
   db_flow_instance = my_houston.db.models.Flow(**flow.to_flow_dict())
   my_houston.db.session.add(db_flow_instance)
   my_houston.db.session.commit()
   print('Has flow key:', db_flow_instance.key.startswith('flow:'))

.. testoutput:: flow_db_test_group

   Has flow key: True

~~~~~~~~~~~~~~
Querying Flows
~~~~~~~~~~~~~~

We can query flows using SqlAlchemy queries via Houston's db utility.

.. testcode:: flow_db_test_group

   flow_from_db = (
      my_houston.db.session.query(my_houston.db.models.Flow)
      .first()
   )
   print('Has flow key:', flow_from_db.key.startswith('flow:'))

.. testoutput:: flow_db_test_group

   Has flow key: True

=====
Tasks
=====
Tasks are an essential component of Flows. See :doc:`tasks`.
