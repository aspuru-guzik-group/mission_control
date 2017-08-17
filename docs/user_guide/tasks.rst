Tasks
=====

===============
What Is a Task?
===============
A MissionControl Task represents a step in a workflow.

A task can do many things. In some cases a task
runs a computation. For example, a task that transforms a JSON string
into a YAML string.

In other cases a task defines logic for how a workflow
should proceed. For example, a task that monitors an
optimization computation and decides whether to keep optimizing or go forward.

In other cases a task wires inputs from one task to another. For example, a task
that that collects the outputs of three separate models and passes them to an
analysis task.

=============
Task Lifecyle
=============
The lifecycle of a task is like this:

#. We define a task as a dict in a flow_spec:

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

#. We convert the flow_spec into a flow and run it with a FlowEngine.
   
#. Eventually the flow reaches our task and the FlowEngine ticks our task:

   .. code-block:: python

     flow_engine.tick_task(task, flow, task_ctx_extras)

#. The FlowEngine calls the tick method of a top-level TaskHandler, and passes
   it the task and the task's context. This context includes the
   task and its parent flow, as well as other user-provided items.

   .. code-block:: python

     task_ctx = {'task': task, 'flow': flow, **task_ctx_extras}
     task_handler.tick_task(task_ctx)

#. The TaskHandler's tick_task method looks up a task-specific tick function
   and calls it with the task_ctx.

#. The tick function performs actions based on the task_ctx.
   The tick method often updates a task's status to indicate whether the task 
   is still running, or whether it has failed or completed.

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
successors is specified when a task is added to a flow, then the task will
have its precursor set to be the task most recently added to the flow. This
default behavior makes it easier to concisely define simple linear flows.

=============
Running Tasks
=============
Now that we can define tasks we can think about how to run them. 

We use tick functions to run tasks. A tick function's only requirement is to
have a signature like this:

.. code-block:: python

   my_tick_fn(*args, task_ctx=None, **kwargs):
       ...

Where 'task_ctx' is a dict that typically has 'task', 'flow', and other
user-specified items.

Per task_type_ above, MissionControl defines a default convention
for loading tick functions based on values in task_ctx.

==============
Tasks Handlers
==============
There are often reasons to define a tick function as a method of a class.

Some of these reasons are:

#. the logic for ticking a task is better expressed as a collection of
   functions, rather than as one long function.

#. multiple tick functions share common setup/teardown logic.

#. multiple tick functions share common logic for handling exceptions.

#. multiple tick functions share common helper utilities and shortcuts.

MissionControl calls a class which defines a tick function a 'TaskHandler'.

MissionControl provides a few classes which you can use as
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
       def initial_tick(self, *args, **kwargs):
           self.task['data']['countdown'] = (
               self.task['task_params']['countdown_start']
           )
           print("Starting countdown for task with key '%s'" % self.task['key'])
           self.intermediate_tick(*args, **kwargs)
       
       def intermediate_tick(self, *args, **kwargs):
           countdown = self.task['data']['countdown']
           if countdown > 0:
               print("T minus %s" % countdown)
               self.task['data']['countdown'] -= 1
           else:
               print("Houston, we have lift off.")
               self.task['data']['my_result'] = 'Lift off!'
               self.task['status'] = 'COMPLETED'

   my_countdown_task = {
       'key': 'my_task',
       'task_params': {
           'countdown_start': 3
       },
       'data': {},
   }

   while my_countdown_task.get('status') != 'COMPLETED':
       CountdownTaskHandler.tick_task(task_ctx={'task': my_countdown_task})

Expected output:

.. testoutput:: task_handler_test_group

   Starting countdown for task with key 'my_task'
   T minus 3
   T minus 2
   T minus 1
   Houston, we have lift off.


=====================
Built-In TaskHandlers
=====================
MissionControl provides a few built-in TaskHandlers for common operations:

.. autoclass:: mc.task_handlers.flow_task_handler.FlowTaskHandler

.. autoclass:: mc.task_handlers.job_task_handler.JobTaskHandler

.. autoclass:: mc.task_handlers.log_task_handler.LogTaskHandler

.. autoclass:: mc.task_handlers.mc_default_task_handler.PrintTaskHandler

.. autoclass:: mc.task_handlers.mc_default_task_handler.NoOpTaskHandler

.. autoclass:: mc.task_handlers.spread_task_handler.SpreadTaskHandler

.. autoclass:: mc.task_handlers.switch_task_handler.SwitchTaskHandler

======
Wiring
======
How do we share data betwen tasks?

What do we mean by wiring? We mean transferring data. But this needs some more
explanation.

For example, say we want to have a flow like this:

#. Get records from a database.

#. For each record: create a subflow.

We could create (1) a task that gets the database records, and (2) a task that
creates subflows. But how would we connect the output of our database task to
our subflow task?

Let's think about some of the strategies we could use.

---------------------------------------
Strategy A: Standard Input/Output Hooks
---------------------------------------
In this strategy we have a convention of using standard input/output hooks.
That is, every
task could have an 'inputs' property and and 'outputs' property. Then, when a
task finishes, a FlowEngine would pass whatever is in the finished task's
'outputs' property to the next task.

This strategy is what many unix command-line tools use. A command takes input
on stdin and then produces output on stdout and stderr.

~~~~~
PROS:
~~~~~
#. Conceptually simple.
#. Works well for simple linear flows.

~~~~~
CONS:
~~~~~
#. Hard to get inputs from multiple sources.
#. Hard to get inputs from previous tasks which are not direct precursors.
#. Obscures how outputs from one task are wired to inputs. Transformations
   must be described within a task's tick function, rather than in a task's
   definition.
#. We need to make intermediate 'adapter' tasks to transform outputs from one
   task into a suitable input format for later tasks.
#. Tasks have to choose what to expose as outputs. If a task does not expose
   a value as an output, a later task cannot access that value.


.. 

------------------------------------------
Strategy B: Tasks Push Inputs To Sucessors
------------------------------------------
In this strategy a task sets inputs on its successors. That is, a task has
references to its successors, and it can set values on its successors.

~~~~~
PROS:
~~~~~
#. No need for adapter tasks. A task controls exactly what gets passed
   to its successors.

~~~~~
CONS:
~~~~~
#. A task has to know what inputs its successors expect. This leads to a high
   degree of coupling, making it hard to reuse or modify tasks.

#. Tasks have to receive or be able to lookup references to their successors.
   This couples tasks to the lookup system.

.. _wiring_strategy:

-------------------------------------
Strategy C: Intermediate Wiring Tasks
-------------------------------------
In this strategy we have two categories of tasks: 'work' tasks that do the work
of our flow, and 'wiring' tasks that pipe data between 'work' tasks.

Our work tasks are ignorant and don't know where their inputs come from.
And they don't make any decisions about what outputs to expose. Their logic only
deals with how to do their work.

On the other hand, our wiring tasks are knowledgeable. They know which
work tasks they pull data from, and what shape of data those source work tasks
expose. Our wiring tasks also know which work tasks they push data to, and what
data formats those target tasks expect.

There is coupling between the wiring tasks and work tasks,
but not between work tasks and other work tasks.  This means that work tasks
can be reused and modified without having to worry about other work tasks.

In addition, our wirings are explicit in the definitions of the wiring tasks.
So we can see them in our flows.

~~~~~
PROS:
~~~~~
#. Work Task Simplicity: A work task's definition can focus exclusively on the
   work its need to do; work tasks no longer need to worry about formatting
   inputs and setting outputs.

#. Decoupling: work tasks can be reused and modified without worrying about
   breaking other work tasks.

#. Multiple Input Sources: a wiring task can aggregate data from multiple
   work tasks into one input.

#. Explicit Connections in Flow Definitions: a wiring task's definition shows
   how data moves through flow, making our flows easier to understand.

~~~~~
CONS:
~~~~~
#. We need to define the intermediate wiring tasks, which can add verbosity.

----------------------
Which Strategy to Use?
----------------------
The MissionControl authors recommend the wiring strategy. We think that this
strategy gives the most control, and makes it easier to reuse and compose
tasks. In addition, we think wiring tasks make flows more comprehensible by
explicitly showing how data moves from one task to another.

To help with creating wiring tasks, MissionControl provides a special
TaskHandler,
the :class:`WireTaskHandler <mc.task_handlers.wire_task_handler.WireTaskHandler>`.

To learn more about wirings see the documentation for  
:class:`WireTaskHandler <mc.task_handlers.wire_task_handler.WireTaskHandler>`, 
and look at the wiring examples.

.. testcode::

  raise NotImplementedError('link to wiring examples')


==============
Interpolations
==============
The primary con of the wiring strategy is that it adds verbosity to a flow.

For example, a simple wiring like 'connect task_a.prop_1 to task_b.prop_2'
requires an entire task definition:

.. testcode::

  raise NotImplementedError('make verbose definition')


MissionControl provides a utility called 'interpolation' to help reduce this
verbosity for simple wirings.

.. testcode::

  raise NotImplementedError('Make inline interpolation example')

The logic for interpolation is defined in the
:class:`mc.task_handlers.mc_default_task_handler.McDefaultTaskHandler` class.
See the :meth:`mc.task_handlers.mc_default_task_handler.McDefaultTaskHandler.get_interpolated_task_ctx` method,
and the interpolation examples.

.. testcode::

  raise NotImplementedError('Link to interpolation examples')


==============
Proxying Tasks
==============
Sometimes we want to create tasks that involve complex definitions. Or we want
to create tasks which change their behavior based on the parameters they
receive.

We could create such tasks by making large, complex task dicts. But these task
dicts can be hard to understand. And if we want to reuse the logic of these
tasks in different flows, we would need to copy-and-paste their task dicts
into our flows. 

Is there a better way? Yes! Proxying tasks to the rescue!

Proxying tasks are tasks which act as simplified interfaces to an underlying
complex task. We can pass simple parameters to a proxying task, and then have
the proxying task dynamically create an underlying complex task from our
parameters.

The way we create a proxying task is by creating a normal task with a special 
dict item named 'proxied_task'. The value of this item is another task dict.

When our FlowEngine encounters a task that has a 'proxied_task' component it
will tick the proxied_task instead.

Let's see an example.

----------------------
Proxying Tasks Example
----------------------

Here's a task with a little complexity. It prints (1) a given message and (2)
the number of characters in the message.

Let's first implement without using a proxying task.

.. testcode:: proxying_tasks_test_group

   from mc.task_handlers.base_task_handler import BaseTaskHandler

   class MessageStatsTaskHandler(BaseTaskHandler):
       def initial_tick(self, *args, **kwargs):
           orig_msg = self.task['task_params']['msg']
           altered_msg = orig_msg + " (%s chars)" % len(orig_msg)

           # Add an additional task to print the altered message.
           flow = self.task_ctx['flow']
           flow.add_task({
               'task_type': 'print',
               'task_params': {
                 'msg': altered_msg
               },
               'precursors': [
                   *self.task['precursors'], self.task['key']
               ]
           })
           self.task['status'] = 'COMPLETED'
   
   flow_spec = {
       'tasks': [
           {
               'task_type': 'message_stats',
               'task_params': {'msg': 'I am short msg'},
           },
           {
               'task_type': 'message_stats',
               'task_params': {'msg': 'I am looooooooong msg'},
           },
       ]
   }

   from mc.flows.flow_engine import FlowEngine
   my_flow_engine = FlowEngine()
   my_flow_engine.run_flow(
       flow=my_flow_engine.flow_spec_to_flow(flow_spec),
       tick_kwargs={
           'tick_fn_for_task_ctx_overrides': {
               'message_stats': MessageStatsTaskHandler.tick_task
            }
       },
   )
    
Expected output:

.. testoutput:: proxying_tasks_test_group

    I am short msg (14 chars)
    I am looooooooong msg (21 chars)

Our implementation feels a bit kludgy. We're altering the structure of our flow
by adding a new task. This is ok for this example. But such an alteration could 
cause problems in other cases. What if other tasks needed to wait until our
altered message was printed? How would they know to use our newly added task as
a successor? How would they be able to get data from it?

Here's a different implementation of our example above. What's different is
that we using a proxying task.

.. testcode:: proxying_tasks_test_group

   class ProxyingMessageStatsTaskHandler(BaseTaskHandler):
       def initial_tick(self, *args, **kwargs):
           orig_msg = self.task['task_params']['msg']
           altered_msg = orig_msg + " (%s chars)" % len(orig_msg)

           self.task['proxied_task'] = {
               'task_type': 'print',
               'task_params': {'msg': altered_msg}
           }

       def intermediate_tick(self, *args, **kwargs):
           self.task['status'] = self.task['proxied_task'].get('status')

   my_flow_engine.run_flow(
       flow=my_flow_engine.flow_spec_to_flow(flow_spec),
       tick_kwargs={
           'tick_fn_for_task_ctx_overrides': {
               'message_stats': ProxyingMessageStatsTaskHandler.tick_task
            }
       },
   )

Expected output:

.. testoutput:: proxying_tasks_test_group

    I am short msg (14 chars)
    I am looooooooong msg (21 chars)

This feels cleaner! The structure of our flow doesn't change, and we don't have
to interact with the flow inside our initial_tick method.

We can trim this down even more if we inherit from MissionControl's
:class:`mc.task_handlers.base_proxying_task_handler.BaseProxyingTaskHandler`
class. This class defines a hook for creating the proxied task, and tasks care
of propagating attributes from the proxied task to the proxying task.

.. testcode: proxying_tasks_test_group

   from mc.task_handlers.base_proxying_task_handler import BaseProxyingTaskHandler

   class SubclassingMessageStatsTaskHandler(BaseProxyingTaskHandler):
       def generate_proxied_task(self):
           orig_msg = self.task['task_params']['msg']
           altered_msg = orig_msg + " (%s chars)" % len(orig_msg)
           return {'task_type': 'print',
                   'task_params': {'msg': altered_msg}}


   my_flow_engine.run_flow(
       flow=my_flow_engine.flow_spec_to_flow(flow_spec),
       tick_kwargs={
           'tick_fn_for_task_ctx_overrides': {
               'message_stats': SubclassingMessageStatsTaskHandler.tick_task
            }
       },
       
Expected output:

.. testoutput:: proxying_tasks_test_group

    I am short msg (14 chars)
    I am looooooooong msg (21 chars)

===============================
Recommended Practices for Tasks
===============================

#. Use TaskHandler classes for non-trivial tick functions.

#. Test your tick functions and task handlers.

#. Use proxying tasks to abstract complex task logic and create reusable tasks.
