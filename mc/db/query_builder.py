import sqlalchemy as _sqla


class QueryBuilder(object):
    CONJUNCTIONS = {
        'AND': {'op': _sqla.and_},
        'OR': {'op': _sqla.or_},
    }

    OPERATORS = {
        '=': lambda f, a: f.__eq__(a),
        '<': lambda f, a: f.__lt__(a),
        '>': lambda f, a: f.__gt__(a),
        '<=': lambda f, a: f.__le__(a),
        '>=': lambda f, a: f.__ge__(a),
        'in': lambda f, a: f.in_(a),
        'ends_with': lambda f, a: f.like('%' + a),
        'starts_with': lambda f, a: f.like(a + '%'),
        'contains': lambda f, a: f.contains(a),
        'is_null': lambda f: f.is_(None),
        'between': lambda f, a: f.between(a[0], a[1])
    }

    def __init__(self, operators=None):
        self.operators = operators or self.OPERATORS

    def alter_query_per_filters(self, query=None, filters=None):
        clause = self._generate_clause_for_filter(
            query=query,
            filter_={
                'conjunction': 'AND',
                'filters': filters
            }
        )
        return query.filter(clause)

    def _generate_clause_for_filter(self, query=None, filter_=None):
        if 'conjunction' in filter_:
            clause = self._generate_clause_for_filter_group(
                query=query, filter_group=filter_)
        else:
            clause = self._generate_clause_for_single_filter(
                query=query, filter_=filter_)
        return clause

    def _generate_clause_for_filter_group(self, query=None,
                                          filter_group=None):
        conjunction = filter_group.get('conjunction')
        if conjunction not in self.CONJUNCTIONS:
            raise Exception("Invalid conjunction '{}'".format(conjunction))
        conjunction_op = self.CONJUNCTIONS[conjunction]['op']
        filter_clauses = [
            self._generate_clause_for_filter(
                query=query,
                filter_=filter_
            )
            for filter_ in filter_group['filters']
        ]
        clause = conjunction_op(*filter_clauses)
        return clause

    def _generate_clause_for_single_filter(self, query=None, filter_=None):
        op = filter_['op']
        negate = False
        if op.startswith('!'):
            negate = True
            op = op.lstrip('! ')
        expr = self._get_expr_for_field(
            query=query, field=filter_['field'])
        args = [expr]
        if 'arg' in filter_:
            args.append(filter_['arg'])
        op_fn = self.operators[op]
        clause = op_fn(*args)
        if negate:
            clause = _sqla.not_(clause)
        return clause

    def _get_expr_for_field(self, query=None, field=None):
        field_components = field.split('.')
        try:
            root_expr = next(col['expr'] for col in query.column_descriptions
                             if col['name'] == field_components[0])
        except StopIteration:
            root_expr = query.column_descriptions[0]['expr']
        if len(field_components) == 1:
            expr = getattr(root_expr, field_components[0])
        else:
            for component in field_components[1:]:
                expr = getattr(expr, component)
        return expr

    def alter_query_per_limit(self, query=None, limit=None):
        return query.limit(limit)

    def alter_query_per_order_by(self, query=None, order_by=None):
        order_by_arg = self._generate_order_by_arg(query=query,
                                                   order_by_spec=order_by)
        return query.order_by(order_by_arg)

    def _generate_order_by_arg(self, query=None, order_by_spec=None):
        expr = self._get_expr_for_field(
            query=query, field=order_by_spec['field'])
        direction = order_by_spec.get('direction') or 'asc'
        order_by_arg = getattr(expr, direction)
        return order_by_arg
