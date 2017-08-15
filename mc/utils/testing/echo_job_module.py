from pathlib import Path
import subprocess
import tempfile
import textwrap


class Constants:
    OUTPUT_FILE_NAME = 'echo_job.out'


def build_work_dir(*args, params=None, output_dir=None, **kwargs):
    output_dir = output_dir or tempfile.mkdtemp()
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    entrypoint_name = 'entrypoint.sh'
    entrypoint_content = textwrap.dedent(
        '''
        #!/bin/bash
        DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
        pushd $DIR > /dev/null
        echo {message} | tee {output_file}
        popd > /dev/null
        '''
    ).lstrip().format(
        message=params['message'],
        output_file=Constants.OUTPUT_FILE_NAME
    )
    entrypoint_path = Path(output_dir, entrypoint_name)
    with open(str(entrypoint_path), 'w') as f:
        f.write(entrypoint_content)
    entrypoint_path.chmod(0x775)
    return {'dir': output_dir, 'entrypoint_name': entrypoint_name}


def run_work_dir(*args, job=None, work_dir=None, use_prebaked=None, **kwargs):
    if use_prebaked:
        with Path(work_dir, Constants.OUTPUT_FILE_NAME).open('w') as f:
            f.write('prebaked')
    else:
        run_cmd = ['bash', str(Path(work_dir, 'entrypoint.sh'))]
        subprocess.run(run_cmd)


def parse_work_dir(*args, job=None, work_dir=None, **kwargs):
    with Path(work_dir, Constants.OUTPUT_FILE_NAME).open() as f:
        output = f.read()
    return {'output': output}
