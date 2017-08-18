EntityDB
========

EntityDB is a framework for storing and querying data in a SQL database.

It provides a SQL implementation of an `EAV Data Model <https://en.wikipedia.org/wiki/Entity%E2%80%93attribute%E2%80%93value_model>`_ . It uses SqlAlchemy to provide python classes and utilities
for creating and querying entity records

=============
Why EntityDB?
=============
MissionControl includes EntityDB for several reasons:

#. Jobs and Flows often need to store and query intermediate data in a variety
   of schemas. EntityDB provides a way to accomodate a wide range of schemas .

#. Jobs and Flows often query data from pre-existing libraries of inputs.
   EntityDB can be used to quickly setup input libraries for testing.

#. It is often useful to decouple data parsing from data ingestion. EntityDB
   provides a system for defining abstract ingestion 'actions', so that job
   parsers can run without needing database connections.

======
Schema
======
The key to understanding how EntityDB works is to understand its schema. That
is, the way it organizes data.

---
Ent
---
The central component of EntityDB's schema is the :class:`mc.db.models.Ent`
class. An Ent represents a 'thing'.

However, when we think about data, we
usually think in more specific terms.  We think about things specific using
terms specific to our field.  For example, if we're thinking
about chemistry, we think in terms of 'molecules', 'conformers',
'electronic structure', and other chemistry-specific terms.

So how do we go from a general 'thing' to a field-specific term? We use Ent's
:attr:`mc.db.models.Ent.ent_type` attribute to indicate that our entity has a
particular type. For example, we could create an Ent that represents a molecule.
In the example below we create a :doc:`houston` instance with an in-memory db,
and create an entity that represents a molecule:

.. testcode :: ent_test_group

  # Setup a Houston instance for easy db access.
  from mc.houston import Houston
  houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
  houston.db.ensure_tables()

  my_molecule = houston.db.models.Ent(ent_type='molecule')
  houston.db.session.add(my_molecule)
  houston.db.session.commit()
  print(my_molecule)

.. testoutput :: ent_test_group
   
   <Ent|ID:...|KEY:ent:molecule:...>

-----
Props
-----
Usually things in our field have different properties. For example, a molecule
can have atoms, or a charge.

We use Ent's :attr:`mc.db.models.Ent.props` attribute to store properties
related to our entity. For example:

.. testcode :: ent_test_group

  my_molecule.props['atoms'] = [
     {'element': 'C', 'x': 0, 'y': 0, 'z': 0},
     {'element': 'O', 'x': 1, 'y': 0, 'z': 0},
     # ...
  ]
  my_molecule.props['charge'] = 3
  houston.db.session.add(my_molecule)
  houston.db.session.commit()
  import json
  print(json.dumps(dict(my_molecule.props), indent=2, sort_keys=True))

.. testoutput :: ent_test_group

   {
     "atoms": [
       {
         "element": "C",
         "x": 0,
         "y": 0,
         "z": 0
       },
       {
         "element": "O",
         "x": 1,
         "y": 0,
         "z": 0
       }
     ],
     "charge": 3
   }

In general you can store anything in props that can be serialized into JSON.
This gives you great flexibility!

.. warning:: 

  If you try to store something in props that can't be serialized you will
  get an error when you try to save your ent.
  
  Example:

  .. testcode :: ent_test_group

     class NonSerializable: pass

     ent_w_bad_props = houston.db.models.Ent(
         props={'this_wont_work': NonSerializable()}
     )
     houston.db.session.add(ent_w_bad_props)
     try:
         houston.db.session.commit()
     except Exception as exc:
         print(exc)

  .. testoutput :: ent_test_group

     (builtins.TypeError) Object of type 'NonSerializable' is not JSON serializable ...

~~~~~~~~~~~~~~~~~~
Filtering by Props
~~~~~~~~~~~~~~~~~~
You can also filter ents by prop values.

Note that we filter on an Ent's :attr:mc.db.models.Ent.props_set attribute, 
rather than :attr:mc.db.models.Ent.props.

Here's an example of filter ents that have
a property with a specific value.

.. testcode :: ent_prop_filter_test_group

  # Setup a Houston instance with a db.
  from mc.houston import Houston
  houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
  houston.db.ensure_tables()

  # shortcut, to save some typing :p
  Ent = houston.db.models.Ent

  # Create some molecule ents.
  molecules = []
  for i in range(1, 4):
      molecule = Ent(
          ent_type='molecule',
          props={'num_atoms': i}
      )
      molecules.append(molecule)
  houston.db.session.add_all(molecules)
  houston.db.session.commit()

  # Setup a base query to filter for molecule ents.
  base_molecule_query = (
       houston.db.session.query(Ent)
      .filter(Ent.ent_type == 'molecule')
  )

  # Filter for molecules with one atom.
  molecules_w_one_atom = (
      base_molecule_query
      .filter(Ent.props_set.any(key='num_atoms', value=1))
      .all()
  )
  for molecule in molecules_w_one_atom:
     print('num_atoms:', molecule.props['num_atoms'])

.. testoutput :: ent_prop_filter_test_group

   num_atoms: 1

Here's a more advanced example that filters for a property matching a numerical
comparison. Notice how we 'join' to Ent.Prop.

.. testcode :: ent_prop_filter_test_group

  # Filter for molecules with greater than one atom.
  molecules_w_multiple_atoms = (
      base_molecule_query
      .join(Ent.Prop, aliased=True, from_joinpoint=True)
      .filter(Ent.Prop.key == 'num_atoms')
      .filter(Ent.Prop.value > 1)
      .order_by(Ent.Prop.value)
      .reset_joinpoint()
      .all()
  )
  for molecule in molecules_w_multiple_atoms:
     print('num_atoms:', molecule.props['num_atoms'])

.. testoutput :: ent_prop_filter_test_group

  num_atoms: 2
  num_atoms: 3

Because EntityDB uses SqlAlchemy, you have great flexibility in the queries
you can build.

See http://docs.sqlalchemy.org/en/latest/orm/query.html#sqlalchemy.orm.query.Query.join
to understand the 'join' line.


-----
Tags
-----

Often we want a way to group ents into collections, or to quickly find
specific ents. EntityDB provides a tagging mechanism to help us do these things.

Here's an example:

.. testcode :: ent_tags_test_group

  # Setup a Houston instance with a db.
  from mc.houston import Houston
  houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
  houston.db.ensure_tables()

  # shortcut, to save some typing :p
  Ent = houston.db.models.Ent

  # Create a molecule ent with some tags.
  # We define tags as python set().
  molecule = Ent(
     ent_type='molecule',
     tags={'some_tag', 'some_other_tag'}
  )
  houston.db.session.add(molecule)
  houston.db.session.commit()
  print(sorted(molecule.tags))

.. testoutput :: ent_tags_test_group

   ['some_other_tag', 'some_tag']

.. testcode :: ent_tags_test_group

   # We can add tags without worrying about duplicates.

   molecule.tags.add('a new tag')
   molecule.tags.add('some_tag')
   houston.db.session.add(molecule)
   houston.db.session.commit()
   print(sorted(molecule.tags))

.. testoutput :: ent_tags_test_group
   
   ['a new tag', 'some_other_tag', 'some_tag']

~~~~~~~~~~~~~~~~~
Filtering by Tags
~~~~~~~~~~~~~~~~~

We can filter ents by tags.

Note that we filter on an Ent's :attr:mc.db.models.Ent.tags_set attribute, 
rather than :attr:mc.db.models.Ent.tagas.

Here's an example of filtering for ents that have
a specific tag.

.. testcode :: ent_filter_tags_test_group

  # Setup a Houston instance with a db.
  from mc.houston import Houston
  houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
  houston.db.ensure_tables()

  # shortcut, to save some typing :p
  Ent = houston.db.models.Ent

  # Create ents with combinations of 2 tags.
  tags = ['tag_%s' % i for i in range(1, 4)]
  molecules = []
  import itertools
  pairs = list(itertools.combinations(tags, 2))
  trios = list(itertools.combinations(tags, 3))
  for tag_combo in [*pairs, *trios]:
      molecule = Ent(
          ent_type='molecule',
          tags=set(tag_combo)
      )
      molecules.append(molecule)
  houston.db.session.add_all(molecules)
  houston.db.session.commit()

  # Filter for ents that have tag_1.
  ents_w_tag_1 = (
     houston.db.session.query(Ent)
     .filter(Ent.tags_set.any(name=tags[0]))
     .all()
  )

  def generate_key_for_tags(tags):
      return '|'.join(sorted(tags))

  tag_set_keys = [generate_key_for_tags(ent.tags) for ent in ents_w_tag_1]
  print(sorted(tag_set_keys))

.. testoutput :: ent_filter_tags_test_group

   ['tag_1|tag_2', 'tag_1|tag_2|tag_3', 'tag_1|tag_3']

In the following example we filter for ents that have a set of tags.

.. testcode :: ent_filter_tags_test_group

  ents_w_tags_1_and_2 = (
     houston.db.session.query(Ent)
     .filter(
         Ent.tags_set.any(name=tags[0])
         & Ent.tags_set.any(name=tags[1])
      )
     .all()
  )

  def generate_key_for_tags(tags):
      return '|'.join(sorted(tags))

  tag_set_keys = [generate_key_for_tags(ent.tags) for ent in ents_w_tags_1_and_2]
  print(sorted(tag_set_keys))

.. testoutput :: ent_filter_tags_test_group

   ['tag_1|tag_2', 'tag_1|tag_2|tag_3']

In the following example we filter for ents that lack a tag.

.. testcode :: ent_filter_tags_test_group

  ents_sans_tag_1 = (
     houston.db.session.query(Ent)
     .filter(~(Ent.tags_set.any(name=tags[0])))
     .all()
  )

  def generate_key_for_tags(tags):
      return '|'.join(sorted(tags))

  tag_set_keys = [generate_key_for_tags(ent.tags) for ent in ents_sans_tag_1]
  print(sorted(tag_set_keys))

.. testoutput :: ent_filter_tags_test_group

   ['tag_2|tag_3']

-------
Lineage
-------
You can organize ents in relationship hierarchies using the attributes

- :attr:`mc.db.models.Ent.parents`
- :attr:`mc.db.models.Ent.children`
- :attr:`mc.db.models.Ent.ancestors`
- :attr:`mc.db.models.Ent.descendants`

Let's see an example.

.. testcode :: ent_lineage_test_group

   # Setup a Houston instance with a db.
   from mc.houston import Houston
   houston = Houston(cfg={'MC_DB_URI': 'sqlite://'})
   houston.db.ensure_tables()

   # shortcut, to save some typing :p
   Ent = houston.db.models.Ent

   # Setup helper methods

   def create_families(num_families=2):
       families = {}
       for i in range(1, 1 + num_families):
           family_key = ('family_%s' % i)
           families[family_key] = create_family(family_key=family_key)
           return families

   def create_family(family_key=None):
       common_props = {'family_key': family_key}
       grandparents = [
           Ent(
               key=('%s:grandparent_%s' % (family_key, i)),
               props={**common_props, 'generation': 'grandparents'}
           )
           for i in range(4)
       ]
       grandparent_pairs = [
           [grandparents[0], grandparents[1]],
           [grandparents[2], grandparents[3]]
       ]
       parents = []
       for i, grandparent_pair in enumerate(grandparent_pairs):
           parents.append(
               Ent(
                   key=('%s:parent_%s' % (family_key, i)),
                   props={**common_props, 'generation': 'parents'},
                   parents=grandparent_pair,
                   ancestors=grandparent_pair
               )
           )
       children = [
           Ent(
               key=('%s:child_%s' % (family_key, i)),
               props={**common_props, 'generation': 'children'},
               parents=parents,
               ancestors=(grandparents + parents)
           )
           for i in range(3)
       ]
       houston.db.session.add_all(grandparents + parents + children)
       houston.db.session.commit()
       family = {
           'grandparents': grandparents,
           'grandparent_pairs': grandparent_pairs,
           'parents': parents,
           'children': children
       }
       return family

   def keys_for_ents(ents): return sorted(set([ent.key for ent in ents]))

   families = create_families(2)
   family_1 = families['family_1']
   individuals = {
       'grandparent': family_1['grandparents'][0],
       'parent': family_1['parents'][0],
       'child': family_1['children'][0]
   }
   for generation, individual in individuals.items():
       for attr in ['parents', 'children', 'ancestors', 'descendants']:
           label = '{generation}.{attr}: '.format(
               generation=generation, attr=attr)
           key_set = keys_for_ents(getattr(individual, attr))
           print(label, key_set)

.. testoutput :: ent_lineage_test_group

   grandparent.parents:  []
   grandparent.children:  ['family_1:parent_0']
   grandparent.ancestors:  []
   grandparent.descendants:  ['family_1:child_0', 'family_1:child_1', 'family_1:child_2', 'family_1:parent_0']
   parent.parents:  ['family_1:grandparent_0', 'family_1:grandparent_1']
   parent.children:  ['family_1:child_0', 'family_1:child_1', 'family_1:child_2']
   parent.ancestors:  ['family_1:grandparent_0', 'family_1:grandparent_1']
   parent.descendants:  ['family_1:child_0', 'family_1:child_1', 'family_1:child_2']
   child.parents:  ['family_1:parent_0', 'family_1:parent_1']
   child.children:  []
   child.ancestors:  ['family_1:grandparent_0', 'family_1:grandparent_1', 'family_1:grandparent_2', 'family_1:grandparent_3', 'family_1:parent_0', 'family_1:parent_1']
   child.descendants:  []

~~~~~~~~~~~~~~~~~~~~
Filtering by Lineage
~~~~~~~~~~~~~~~~~~~~
We can filter ents by their lineage.

Note that unlike props or tags, we filter on the direct lineage attributes. The
reason for this is technical: props and tags are SqlAlchemy association_proxies,
whereas lineage attributes are SqlAlchemy relationships.

Here are some examples.

Querying by parents:

.. testcode :: ent_lineage_test_group

   grandparent_pair_0 = family_1['grandparent_pairs'][0]
   children_of_grandparent_pair_0 = (
       houston.db.session.query(Ent)
       .join(Ent.parents, aliased=True, from_joinpoint=True)
       .filter(
           Ent.key.in_([
               grandparent.key
               for grandparent in grandparent_pair_0
           ])
       )
       .reset_joinpoint()
       .all()
   )
   print("\n".join(keys_for_ents(children_of_grandparent_pair_0)))

.. testoutput :: ent_lineage_test_group

   family_1:parent_0

Querying by ancestors:

.. testcode :: ent_lineage_test_group

   descendants = (
       houston.db.session.query(Ent)
       .filter(
           Ent.props_set.any(key='generation', value='children')
       )
       .join(Ent.ancestors, aliased=True, from_joinpoint=True)
       .filter(
           Ent.props_set.any(key='family_key', value='family_1')
       )
       .reset_joinpoint()
       .all()
   )
   print("\n".join(keys_for_ents(descendants)))

.. testoutput :: ent_lineage_test_group

    family_1:child_0
    family_1:child_1
    family_1:child_2


======================================
Updating an EntityDB via Actions
======================================
EntityDB allows you to describe updates to a db as a list of dicts. We call
these dicts 'actions'.

An action is a dict with keys for 'type' and 'params'.

Actions are useful because they allow a program to generate a set of updates
without needing a DB connection. This is allows us to do things like:

#. Define parsers for job artifacts that can run without a DB connection.

#. Test parse results by in terms of expected actions, rather than expected DB
   state.

#. Batch together DB writes into chunks, for more efficient write performance.

--------------
Upsert Actions
--------------

Currently the only valid type is 'upsert'.

The general idea of an upsert action is to get or create an ent, and then
run a series of updates to that ent.

See :meth:`mc.db.Db.upsert` for how
the parameters that an upsert action can take.

.. testcode ::

  raise NotImplementedError('link upsert action example')

========================================
Recommended Practices for Using EntityDB 
========================================

#. Think carefully about your schema before you start creating large numbers of
   entities.

   #. Can you design queries that get the correct answers for your questions?
      
   #. What kind of lineage relationships do you need to track?

   #. What kind of properties do you expect to store? Will you need to filter
      on these properties?

#. Test your schema on a small scale before you start creating large numbers of
   entities.
