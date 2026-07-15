import json
import os
from threading import Lock

from output import logger

METADATA_FILE = '.metadata'

_metadata_lock = Lock()


def metadata_path(directory):
    return os.path.join(directory, METADATA_FILE)


def load_metadata(directory):
    try:
        with open(metadata_path(directory)) as f:
            return json.load(f)
    except FileNotFoundError:
        return _legacy_metadata(directory)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Could not read {metadata_path(directory)}: {e}")
        return {}


def update_metadata(directory, **fields):
    with _metadata_lock:
        data = load_metadata(directory)
        data.update(fields)
        path = metadata_path(directory)
        tmp_path = path + '.tmp'
        try:
            with open(tmp_path, 'w') as f:
                json.dump(data, f, indent=2, sort_keys=True)
                f.write('\n')
            os.replace(tmp_path, path)
        except OSError as e:
            logger.warning(f"Could not write {path}: {e}")
    return data


def _legacy_metadata(directory):
    data = {}

    def _first_existing(name):
        for path in (os.path.join(directory, name),
                     os.path.join(directory, 'log_files', name)):
            if os.path.exists(path):
                return path
        return None

    path = _first_existing('.method')
    if path:
        try:
            with open(path) as f:
                data['method'] = f.read().strip()
        except OSError:
            pass

    path = _first_existing('.constrain')
    if path:
        try:
            with open(path) as f:
                indexes = {}
                for line in f:
                    atom, index = line.split(':')
                    indexes[atom.strip()] = int(index.strip())
            # Migrate the pre-Cl/NO3 key names O/OH to X/XH
            if 'O' in indexes and 'X' not in indexes:
                indexes['X'] = indexes.pop('O')
                if 'OH' in indexes:
                    indexes['XH'] = indexes.pop('OH')
            data['constrained_indexes'] = indexes
        except (OSError, ValueError):
            logger.warning(f"Could not parse legacy constraint file {path}")

    path = _first_existing('.symmetry')
    if path:
        try:
            with open(path) as f:
                data['reaction_path_degeneracy'] = int(f.read().strip())
        except (OSError, ValueError):
            pass

    return data
