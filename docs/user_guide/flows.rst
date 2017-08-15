Flows
=====

==============
What Is a Flow?
==============
A MissionControl flow represents a set of tasks run in a specific sequence,
based on the state of other tasks.

============
Flow Lifecyle
============
In general the lifecycle of a flow is like this:

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

#. Tick the flow until it fails, or until all of its tasks are finished:
   ::

     while flow.status not in {'FAILED', 'COMPLETED'}:
         tick_flow(flow)

MissionControl defines conventions and commands to help manage the
flow lifecycle.

==============
Defining Flows
==============
MissionControl defines abstract flow specs as dictionaries with a few primary
components:

label
  a human-readable label for the flow

key
  Often you will want to have a key that uniquely identifies a flow.
  For example, a UUID. This key is useful for tracking flow statuses, and for
  avoiding name collisions when you work with multiple flows.

cfg
  Configuration parameters for how the flow should run. For example, if one
  task fails should the entire flow fail?

data
  Initial data that the flow should have.

tasks
  A list of task_specs. More on task_specs below.

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
Now we can define flow specs. But how do turn our flow_spec into a Flow object?

We use :mod:`FlowEngine<mc.flows.flow_engine.FlowEngine>` . FlowEngine is
a class which contains methods for converting between flow formats and for
ticking flows.

To convert our flow_spec into a flow object we call
:meth:`mc.flows.flow_engine.FlowEngine.flow_spec_to_flow`.

We get back an instance of :class:`mc.flows.flow.Flow`. This class contains
methods for querying flow tasks and manipulating a flow's underlying
attributes.

=============
Running Flows
=============
Now that we have a Flow object, we can run our flow with the FlowEngine.

.. testcode:

    print()

.. testoutput

   foo

==============
Storing Flows
==============

=====
Tasks
=====

============================================
Recommended Practices for Working with Flows
============================================
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
