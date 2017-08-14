from pathlib import Path
import tempfile
import textwrap


def build_work_dir(*args, params=None, output_dir=None, **kwargs):
    output_dir = output_dir or tempfile.mkdtemp()
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    entrypoint_name = 'entrypoint.sh'
    entrypoint_content = textwrap.dedent(
        '''
        #!/bin/bash
        echo {message}
        '''
    ).lstrip().format(message=params['message'])
    entrypoint_path = Path(output_dir, entrypoint_name)
    with open(str(entrypoint_path), 'w') as f:
        f.write(entrypoint_content)
    entrypoint_path.chmod(0x775)
    return {'dir': output_dir, 'entrypoint_name': entrypoint_name}
