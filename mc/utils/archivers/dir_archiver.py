from pathlib import Path
import shutil
import datetime

from .base_archiver import BaseArchiver


class DirArchiver(BaseArchiver):
    PATH_COMPONENT_SEPARATOR = '__'

    def __init__(self, *args, root_dir=None, ensure_root=True,
                 path_components_generator=None,
                 time_components_generator=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.root_path = Path(root_dir)
        self.path_components_generator = (
            path_components_generator or
            self._default_path_components_generator
        )
        self.time_components_generator = (
            time_components_generator or
            self._default_time_components_generator
        )
        if ensure_root:
            self.root_path.mkdir(parents=True, exist_ok=True)

    def _default_path_components_generator(self, src=None):
        path_components = []
        if self.time_components_generator:
            path_components.extend(self.time_components_generator(src=src))
        basename = str(Path(src).name)
        path_components.extend(basename.split(self.PATH_COMPONENT_SEPARATOR))
        return path_components

    def _default_time_components_generator(self, src=None):
        return [datetime.datetime.utcnow().isoformat(timespec='hours')]

    def ingest(self, src=None, transfer_fn=shutil.move):
        rel_dest_path = self._generate_rel_dest_path(src=src)
        abs_dest_path = Path(self.root_path, rel_dest_path)
        transfer_fn(src, str(abs_dest_path))
        meta = {'key': str(rel_dest_path)}
        return meta

    def _generate_rel_dest_path(self, src=None):
        return Path(*self.path_components_generator(src=src))

    def materialize_as_path(self, meta=None):
        return (self.root_path / meta['key'])
