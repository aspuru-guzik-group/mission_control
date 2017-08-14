from pathlib import Path
import tempfile
import textwrap


def build_work_dir(*args, params=None, output_dir=None, **kwargs):
    output_dir = output_dir or tempfile.mkdtemp()
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    entrypoint_name = 'entrypoint.sh'
    with open(Path(output_dir, entrypoint_name), 'w') as f:
        f.write(textwrap.dedent(
            '''
            #!/bin/bash
            echo "noop"
            '''
        ))
    return {'dir': output_dir, 'entrypoint_name': entrypoint_name}


def parse_work_dir(*args, **kwargs): return []
