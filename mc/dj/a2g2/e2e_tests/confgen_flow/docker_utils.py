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
        self.ody_ip_addr = self.get_container_ip_addr(
            network=self.mc_network,
            container_info=container_infos_by_service['odyssey']
        )

    def start_dockers(self):
        subprocess.run(
            'cd {docker_dir} && \
            MC_ROOT={mc_root} \
            MC_NETWORK={mc_network} \
            docker-compose up -d'.format(
                docker_dir=self.docker_dir,
                mc_root=os.environ.get('MC_ROOT', mc.__path__),
                mc_network=self.mc_network,
            ),
            shell=True
        )

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
        container_ids = subprocess.check_output(
            'cd {docker_dir} && \
            MC_NETWORK={mc_network} \
            docker-compose ps -q'.format(
                docker_dir=self.docker_dir,
                mc_network=self.mc_network,
            ),
            shell=True
        ).decode()
        return container_ids.split()

    def key_container_infos_by_service(self, container_infos=None):
        keyed_container_infos = {}
        for container_info in container_infos.values():
            service = container_info['Config']['Labels']\
                    ['com.docker.compose.service']
            keyed_container_infos[service] = container_info
        return keyed_container_infos

    def get_container_ip_addr(self, network=None, container_info=None):
        #print(json.dumps(container_info, indent=2))
        return container_info['NetworkSettings']['Networks']\
                [network]['IPAddress']

    def teardown(self):
        subprocess.run(
            'cd {docker_dir} && \
            docker-compose down'.format(
                docker_dir=self.docker_dir
            ),
            shell=True
        )
