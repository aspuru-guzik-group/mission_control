import os
import string
import random
import subprocess
import pexpect


class NotAuthorizedException(Exception):
    pass

class SSHControlSocketClient(object):
    def __init__(self, user=None, host=None, control_socket=None,
                 control_dir=None):
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

    def _ensure_dir(self, _dir):
        if not os.path.isdir(_dir):
            os.makedirs(_dir)
        return _dir

    def connect(self):
        if not self._is_authorized():
            child = pexpect.spawn("ssh -fNMS %s %s@%s" % (
                self.control_socket, self.user, self.host))
            child.interact()

    def _is_authorized(self, check=False):
        completed_process = self._run_ssh_control_cmd('check', check=check)
        return (completed_process.returncode == 0)

    def disconnect(self):
        if self._is_authorized():
            self._run_ssh_control_cmd('exit')

    def _run_ssh_control_cmd(self, control_cmd, check=False):
        return subprocess.run(
            ['ssh', '-S', self.control_socket, '-O', control_cmd, 'go'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=check)

    def run(self, command, check=False):
        self._ensure_authorized()
        return subprocess.run(
            ['ssh', '-S', self.control_socket, 'go'] + command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=check,
            universal_newlines=True)

    def scp(self, src=None, dst=None, flags=''):
        self._ensure_authorized()
        prefixed_dst = '%s@%s:%s' % (self.user, self.host, dst)
        subprocess.run(['scp', '-o ControlPath=%s' % self.control_socket,
                        flags, src, prefixed_dst], check=True)

    def _ensure_authorized(self):
        self._is_authorized(check=True)
