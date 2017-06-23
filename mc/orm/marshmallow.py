import marshmallow
import sqlalchemy.types as _sa_types

from mc import constants as _mc_constants
from . import custom_sa_types as _custom_sa_types

def sa_schema_to_marsh_schemas(sa_schema=None, Schema=marshmallow.Schema,
                               fields=marshmallow.fields):
    marshmallow_schemas = {}

    def sa_column_to_marshmallow_field(sa_column=None):
        field_cls = None
        field_kwargs = {}
        if type(sa_column.type) is _sa_types.String:
            field_cls = fields.String
            field_kwargs['max_length'] = sa_column.type.length
        elif type(sa_column.type) is _sa_types.Integer:
            field_cls = fields.Integer
        elif type(sa_column.type) is _sa_types.Text:
            field_cls = fields.String
        elif type(sa_column.type) is _sa_types.Boolean:
            field_cls = fields.Boolean
        elif type(sa_column.type) is _sa_types.DateTime:
            field_cls = fields.DateTime
            field_kwargs['format'] = _mc_constants.DATETIME_FORMAT
        elif type(sa_column.type) is _custom_sa_types.JSONEncodedDict:
            field_cls = fields.Dict
        assert field_cls is not None, "no field class set"
        return field_cls(**field_kwargs)

    common_dump_only_fields = ['created', 'modified']

    flow_table = sa_schema['tables']['Flow']
    class FlowSchema(Schema):
        key = sa_column_to_marshmallow_field(flow_table.columns['key'])
        label = sa_column_to_marshmallow_field(flow_table.columns['label'])
        serialization = sa_column_to_marshmallow_field(
            flow_table.columns['serialization'])
        status = sa_column_to_marshmallow_field(flow_table.columns['status'])
        num_tickable_tasks = sa_column_to_marshmallow_field(
            flow_table.columns['num_tickable_tasks'])
        depth = sa_column_to_marshmallow_field(flow_table.columns['depth'])
        claimed = sa_column_to_marshmallow_field(flow_table.columns['claimed'])
        created = sa_column_to_marshmallow_field(flow_table.columns['created'])
        modified = sa_column_to_marshmallow_field(flow_table.columns['modified'])
        class Meta:
            fields = ['key', 'label', 'serialization', 'status', 'claimed',
                      'created', 'modified']
            dump_only = common_dump_only_fields
    marshmallow_schemas['Flow'] = FlowSchema()

    job_table = sa_schema['tables']['Job']
    class JobSchema(Schema):
        key = sa_column_to_marshmallow_field(job_table.columns['key'])
        label = sa_column_to_marshmallow_field(job_table.columns['label'])
        job_spec = sa_column_to_marshmallow_field(job_table.columns['job_spec'])
        data = sa_column_to_marshmallow_field(job_table.columns['data'])
        status = sa_column_to_marshmallow_field(job_table.columns['status'])
        claimed = sa_column_to_marshmallow_field(job_table.columns['claimed'])
        created = sa_column_to_marshmallow_field(job_table.columns['created'])
        modified = sa_column_to_marshmallow_field(job_table.columns['modified'])
        class Meta:
            fields = ['key', 'label', 'job_spec', 'data', 'status', 'claimed',
                      'created', 'modified']
            dump_only_fields = common_dump_only_fields
    marshmallow_schemas['Job'] = JobSchema()

    queue_table = sa_schema['tables']['Queue']
    class QueueSchema(Schema):
        key = sa_column_to_marshmallow_field(queue_table.columns['key'])
        label = sa_column_to_marshmallow_field(queue_table.columns['label'])
        queue_spec = sa_column_to_marshmallow_field(
            queue_table.columns['queue_spec'])
        created = sa_column_to_marshmallow_field(queue_table.columns['created'])
        modified = sa_column_to_marshmallow_field(queue_table.columns['modified'])
        class Meta:
            fields = ['key', 'label', 'queue_spec', 'created', 'modified']
            dump_only_fields = common_dump_only_fields
    marshmallow_schemas['Queue'] = QueueSchema()

    return marshmallow_schemas
