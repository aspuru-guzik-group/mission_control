from flow_engines.flow import Flow

from ..nodes.base import BaseNode
from ..nodes.job import Job_Node


class RPI_Flow_Module(object):
    flow_type = 'rpi'
    label = 'RPI'
    description = 'Flow that runs run-parse-ingest cycle.'

    @classmethod
    def generate_flow(cls, *arg, **kwargs):
        node_classes = cls.get_node_classes()
        flow = Flow()

        run_node = node_classes['rpi:job'](id='run')
        flow.add_node(node=run_node, as_root=True)

        parse_setup_node = node_classes['rpi:parse_setup'](id='parse_setup')
        flow.add_node(node=parse_setup_node, precursor=run_node)
        parse_node = node_classes['rpi:job'](id='parse')
        flow.add_node(node=parse_node, precursor=parse_setup_node)

        ingest_setup_node = node_classes['rpi:ingest_setup'](id='ingest_setup')
        flow.add_node(node=ingest_setup_node, precursor=parse_node)
        ingest_node = node_classes['rpi:job'](id='ingest')
        flow.add_node(node=ingest_node, precursor=ingest_setup_node)


    @classmethod
    def get_node_classes(cls):
        node_classes = {
            'rpi:job': Job_Node,
            'rpi:parse_setup': Parse_Setup_Node,
            'rpi:ingest_setup': Ingest_Setup_Node,
        }
        return node_classes

class Parse_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.flow.nodes['query']
        dest_node = self.flow.nodes['run']
        dest_node.data['input'] = {
            'smiles': src_node.data['output']['result_set']['mols'][0]['smiles'],
            'confgen_spec': '!!!@TODO!!!'
        }
        self.status = 'COMPLETED'

class Ingest_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.flow.nodes['query']
        dest_node = self.flow.nodes['run']
        dest_node.data['input'] = {
            'smiles': src_node.data['output']['result_set']['mols'][0]['smiles'],
            'confgen_spec': '!!!@TODO!!!'
        }
        self.status = 'COMPLETED'
