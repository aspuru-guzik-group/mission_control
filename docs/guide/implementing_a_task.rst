Implementing a Task
-------------------
Now we know we want to have a certain task in our flow. 

How do we implement it?

Again, the key is to break it down without getting bogged down.

A few heuristic questions:

#. Can we break the task down into subtasks?
#. Does our task need access to special **resources**? For example, does it need
#  to read a file or access a database? If yes, it may make sense to have our
#  task create a job.
#. Does our task need access to special **resources**? For example, does it need
#  to read a file or access a database? If yes, it may make sense to have our
#  task create a job.
#. What **inputs** does our task need?
#. What **outputs** should our task produce?
#. Where do we want the code for our task to live?
#. How can we test that our task does what we expect?

