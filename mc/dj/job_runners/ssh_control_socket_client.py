import os
import logging
import string
import random
import subprocess
import pexpect


class NotAuthorizedException(Exception):
    pass

class SSHControlSocketClient(object):
    def __init__(self, user=None, host=None, control_socket=None,
                 control_dir=None, auth_fn=None, logger=None):
        self.user = user or os.environ.get('USER')
        self.host = host
        if not control_socket:
            if not control_dir:
                control_dir = os.path.join(os.environ['HOME'],
                                           '.ssh_control_sockets')
            control_dir = self._ensure_dir(control_dir)
            socket_id = ''.join([random.choice(string.ascii_letters)
                                 for i in range(10)])
            control_socket = os.path.join(control_dir, socket_id)
        self.control_socket = control_socket
        self.logger = logger or logging
        self.auth_fn = auth_fn or self.default_auth_fn

    def _ensure_dir(self, _dir):
        if not os.path.isdir(_dir):
            os.makedirs(_dir)
        return _dir

    def connect(self):
        if not self._is_authorized():
            ssh_cmd = ("ssh -o \"StrictHostKeyChecking no\" -fNMS {socket} "
                       "{user}@{host}"
                      ).format(socket=self.control_socket,
                               user=self.user,
                               host=self.host)
            self.auth_fn(ssh_client=self, ssh_cmd=ssh_cmd)
        self.logger.info("connected")

    def default_auth_fn(self, ssh_cmd=None, **kwargs):
        pexpect.spawn(ssh_cmd).interact()

    def _is_authorized(self, check=False):
        completed_process = self._run_ssh_control_cmd('check', check=check)
        return (completed_process.returncode == 0)

    def disconnect(self):
        if self._is_authorized():
            self._run_ssh_control_cmd('exit')

    def _run_ssh_control_cmd(self, control_cmd, check=False):
        p = subprocess.run(
            ['ssh', '-S', self.control_socket, '-O', control_cmd, 'go'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=check)
        return p

    def run_process(self, cmd=None, check=False):
        self._ensure_authorized()
        return subprocess.run(
            ['ssh', '-S', self.control_socket, 'go'] + cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=check,
            universal_newlines=True)

    def scp(self, src=None, dest=None, flags=''):
        self._ensure_authorized()
        subprocess.run(['scp', '-o ControlPath=%s' % self.control_socket,
                        flags, src, dest], check=True)

    def scp_from_remote(self, remote_src_path=None, local_dest_path=None,
                        flags=''):
        prefixed_src_path = '%s@%s:%s' % (self.user, self.host, remote_src_path)
        return self.scp(src=prefixed_src_path, dest=local_dest_path,
                        flags=flags)

    def scp_to_remote(self, local_src_path=None, remote_dest_path=None,
                      flags=''):
        prefixed_dest_path = '%s@%s:%s' % (self.user, self.host,
                                           remote_dest_path)
        return self.scp(src=local_src_path, dest=prefixed_dest_path,
                        flags=flags)

    def rsync(self, src=None, dest=None, flags=''):
        self._ensure_authorized()
        cmd = ['rsync', '-e', 'ssh -o "ControlPath=%s"' % self.control_socket,
               flags, src, dest]
        subprocess.run(cmd, check=True)

    def rsync_from_remote(self, remote_src_path=None, local_dest_path=None,
                          flags=''):
        prefixed_src_path = '%s@%s:%s' % (self.user, self.host, remote_src_path)
        return self.rsync(src=prefixed_src_path, dest=local_dest_path,
                        flags=flags)

    def rsync_to_remote(self, local_src_path=None, remote_dest_path=None, flags=''):
        prefixed_dest_path = '%s@%s:%s' % (self.user, self.host,
                                           remote_dest_path)
        return self.rsync(src=local_src_path, dest=prefixed_dest_path,
                          flags=flags)

    def _ensure_authorized(self):
        try:
            self._is_authorized(check=True)
        except Exception as error:
            raise Exception("SSH connection not authorized: %s" % error)
