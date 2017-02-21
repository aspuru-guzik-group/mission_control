import json
import os
import subprocess

import mc


class DockerEnv(object):
    def __init__(self):
        this_dir = os.path.dirname(__file__)
        self.docker_dir = os.path.join(this_dir, 'docker')

    def setup(self):
        self.mc_network = self.get_mc_network()
        self.start_dockers()
        container_infos = self.get_container_infos()
        container_infos_by_service = self.key_container_infos_by_service(
            container_infos=container_infos)
        self.mc_ip_addr = self.get_container_ip_addr(
            network=self.mc_network,
            container_info=container_infos_by_service['mission_control']
        )
        self.odyssey_user_username = 'test_user'
        self.odyssey_user_password = 'test_pass'
        self.odyssey_ip_addr = self.get_container_ip_addr(
            network=self.mc_network,
            container_info=container_infos_by_service['odyssey']
        )

    def start_dockers(self):
        return self.run_compose_cmd('up -d')

    def run_compose_cmd(self, compose_cmd=None, check=True, **run_kwargs):
        completed_proc = subprocess.run(
            'cd {docker_dir} && \
            MC_ROOT={mc_root} \
            MC_NETWORK={mc_network} \
            docker-compose {compose_cmd}'.format(
                docker_dir=self.docker_dir,
                mc_root=os.environ.get('MC_ROOT', mc.__path__),
                mc_network=self.mc_network,
                compose_cmd=compose_cmd
            ),
            stdout=subprocess.PIPE,
            shell=True,
            check=check,
            **run_kwargs
        )
        return completed_proc

    def get_mc_network(self):
        mc_network = subprocess.check_output(
            'docker network ls \
            --filter "label=mc.network=mc" \
            --format {{.Name}}',
            shell=True,
        ).decode().strip()
        return mc_network

    def get_container_infos(self):
        container_infos = {
            container_id: self.get_inspection_info(object_id=container_id)
            for container_id in self.get_container_ids()
        }
        return container_infos

    def get_inspection_info(self, object_id=None):
        inspection_output = subprocess.check_output(
            'docker inspect {object_id}'.format(object_id=object_id),
            shell=True
        ).decode()
        inspection_info = json.loads(inspection_output)[0]
        return inspection_info

    def get_container_ids(self):
        container_ids = self.run_compose_cmd('ps -q', check=True)\
                .stdout.decode()
        return container_ids.split()

    def key_container_infos_by_service(self, container_infos=None):
        keyed_container_infos = {}
        for container_info in container_infos.values():
            service = container_info['Config']['Labels']\
                    ['com.docker.compose.service']
            keyed_container_infos[service] = container_info
        return keyed_container_infos

    def get_container_ip_addr(self, network=None, container_info=None):
        return container_info['NetworkSettings']['Networks']\
                [network]['IPAddress']

    def teardown(self):
        self.run_compose_cmd('down', check=True)
