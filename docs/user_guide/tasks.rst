Tasks
=====

==============
What Is a Task?
==============
A MissionControl Task represents a step in a workflow.

A task can do many things. In some cases a task
runs a computation. For example, a task that transforms a JSON string
into a YAML string.

In other cases a task defines logic for how a workflow
should proceed. For example, a task that 
optimization computation and decides whether to keep optimizing or go forward.

In other cases a task wires inputs from one task to another. For example, a task
that that collects the outputs of three separate models and passes them to an
analysis task.

============
Task Lifecyle
============
The lifecycle of a task is like this:

#. We define a task as a dict within in a flow_spec:
   .. code-block:: python

     flow_spec = {
         'label': 'example_flow',
         'tasks': [
             {
                 'key': 'task_1',
                 'task_type': 'print',
                 'task_params': {'msg': 'I am task_1.'},
             },
             ...
         ]
     }

#. We build a flow from the flow_spec and run it with a FlowEngine. Eventually
   the flow reaches a state when it's time for our task to run. At this point
   the FlowEngine ticks our task.
   .. code-block:: python

     flow_engine.tick_task(task, flow, task_ctx)

#. When the FlowEngine ticks our task it calls a top-level TaskHandler. The
   TaskHandler calls its tick method with a 'task_ctx' dict that contains the
   task and its parent flow, as well as other user-provided items.
   .. code-block:: python

     task_handler.tick_task(task_ctx={'task': task, 'flow': flow, ...})

#. The TaskHandler's tick_task method looks up a task-specific tick function
   and calls it.

#. The the tick function performs some set of actions based on the task_ctx.
   The tick method often also updates a task's status to indicate whether it
   is still running, or whether it has failed or completed.

MissionControl defines conventions and provides utilities to help define tasks
and run them.

==============
Defining Tasks
==============
A task is defined as a dict, usually in the context of a flow spec. A task has
several components: 

label
  a human-readable label for the task

key
  A key that uniquely identifies a task. This key is usually human-readable,
  and can be used later to access attributes of the task.
  
  Keys should be formatted as valid python variable names
  (no dots or special characters).
  'my_task' is a valid key, while 'my.task' is not.

task_type
  .. _task_type:

  An identifier for the type of task. MissionControl's default task handler
  will look up an appropriate handler for the task by dispatching on the value
  of task_type.

  From :class:`mc.task_handlers.mc_default_task_handler.McDefaultTaskHandler`:

  .. automethod:: mc.task_handlers.mc_default_task_handler.McDefaultTaskHandler.default_task_type_to_tick_fn_dot_spec

task_params
  A dict of parameters for the task.

data
  A dict data that the task should have.

precursors
  A list of strings referring to keys of tasks which must complete before this
  task can run.

successors
  A list of strings referring to keys of tasks which should not start before
  this task finishes.

Example:

.. code-block:: python

  my_task = {
     'key': 'my_task',
     'task_type': 'my.task.type',
     'task_params': {'some_param': 'some_value'},
     'data': {'some_key': 'some_value'},
     'precursors': ['some_task_that_should_finish_before_my_task'],
     'successors': ['some_task_that_needs_to_wait_for_my_task'],
  }

A special note on precursors and successors: if neither precursors or
successors is specified, then when a task is added to a flow, it will have the
last task added to the flow set as its precursor. This makes it easy to define
simple linear flows of tasks concisely.

=============
Running Tasks
=============
Now that we can define tasks we can think about how to run them. 

We use tick functions to run tasks. The only requirements for a tick function
is to have a signature like this:

.. code-block:: python

   my_tick_fn(*args, task_ctx=None, **kwargs):
       ...

Where 'task_ctx' is a dict that typically has 'task', 'flow', and other
user-specified items.

Per task_type_ above, MissionControl defines a default convention
for loading tick functions.

==============
Tasks Handlers
==============
There are often reasons to define tick functions in the context of a class.
MissionControl calls these classes 'TaskHandlers'.

Some common reasons to define a TaskHandler class instead of a plain tick
function include:

#. the logic for ticking a task is better expressed as a collection of
   functions, rather than as one long function.

#. multiple tick functions share common setup/teardown logic.

#. multiple tick functions share common logic for handling exceptions.

#. multiple tick functions share common helper utilities.

MissionControl provides a few base TaskHandler classes which you can use as
base classes for your own TaskHandlers.

See :class:`mc.task_handlers.base_task_handler` to find out more about these 
base TaskHandler classes.

Usually you will want to inherit from
:class:`mc.task_handlers.base_task_handler.BaseTaskHandler` and then define
implementations for the methods
:meth:`mc.task_handlers.base_task_handler.BaseTaskHandler.initial_tick` and
:meth:`mc.task_handlers.base_task_handler.BaseTaskHandler.intermediate_tick`.

Example:

.. testcode:: task_handler_test_group

   from mc.task_handlers.base_task_handler import BaseTaskHandler

   class CountdownTaskHandler(BaseTaskHandler):
       def initial_tick(self, *args, task_ctx=None, **kwargs):
           self.task['data']['countdown'] = (
               self.task['task_params']['countdown_start']
           )
           print("Starting countdown for task with key '%s'" % self.task['key'])
           self.intermediate_tick(*args, task_ctx=task_ctx, **kwargs)
       
       def intermediate_tick(self, *args, task_ctx=None, **kwargs):
           countdown = self.task['data']['countdown']
           if countdown > 0:
               print("T minus %s" % countdown)
               self.task['data']['countdown'] -= 1
           else:
               print("Houston, we have lift off.")
               self.task['data']['my_result'] = 'Lift off!'
               self.task['status'] = 'COMPLETED'

   my_countdown_task = {
       'task_params': {
           'countdown_start': 3
       },
       'data': {},
   }

   while my_countdown_task.get('status') != 'COMPLETED':
       CountdownTaskHandler.tick_task(task_ctx={'task': my_countdown_task})

Expected output:

.. testoutput:: task_handler_test_group

   T minus 3
   T minus 2
   T minus 1
   Houston, we have lift off.
