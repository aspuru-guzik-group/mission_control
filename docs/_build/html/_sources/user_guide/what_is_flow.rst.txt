What Is a Flow?
===============

A flow is a collection of *tasks* which should execute in a specific *sequence*.

Modeling Flows: Flows as DAGs
-----------------------------

In MissionControl, we think about flows as *Directed Acyclic Graphs*, or *DAGs*.
(For more on DAGs, see <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_.) In our DAG the nodes are our tasks; the edges define the sequence.

Modeling flows as DAG allows us to define many types of flows. For example, we can define a flow that is a simpler linear sequence:

.. todo::

  - MAKE LINEAR FLOW GRAPHIC

Or, we can define a flow that has branches that run in parallel:

.. todo::

  - MAKE PARALEL FLOW GRAPHIC

Representing Flows
------------------
So now we know how we can model flows. But how can we represent a flow? That is, how can we write a flow?

One way to represent a flow is as a drawing. For example:

.. todo::

  - MAKE FLOW DRAWING

Drawings are often the representation we use when we first design a flow. 

Ideally we could make a drawing, and then show that drawing to a computer and have it run our flow.
But one of the problems with drawings is that they are hard for computers to understand.

.. note::

  A challenging exercise for the enterprising reader: make a tool that converts
  a drawing of a flow into a text representation.

Since we want to make our flows readable by computers, we use other representations.

One way we can represent a flow is as a list of tasks and edges.

Another way we can represent a flow is as a list of tasks and edges.

Going From a Drawing to a Flow
------------------------------
How do we go from a drawing to a running flow? We break it down.

There are a few ways to do this breakdown. One way is to make a complete
skeleton of our flow, and then fill it in. Another is way is to add tasks one at
a time, testing them and building them up along the way.

Both approaches have merits, and you can even use a combination of the two.

The key thing is not to get bogged down! Remember, you can always change your
flow later. So just write as much as you know now, and fill in the rest later.

One way is to start by making a skeleton flow that represents the sequencing
of our flow. In this skeleton flow our tasks don't do anything yet. They just print out a message.

.. todo:

   - MAKE SKELETON FLOW EXAMPLE
   - MAKE INCREMENTAL EXAMPLE

This skeleton gives us a starting point. We can check that the sequence of our flow works as expected. And then we can fill logic in as needed.

The next step is implementing a task.

.. todo:

   - MAKE LINK FOR IMPLEMENTING A TASK
