import json
import os
import sys
from threading import Lock

import runtime
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


def guard_monitor_pid():
    if runtime.args.restart or runtime.args.rerun:
        old_pid = load_metadata(runtime.start_dir).get('monitor_pid')
        if old_pid:
            try:
                os.kill(int(old_pid), 0)
            except (OSError, ValueError):
                pass  # unreadable pid or the process is gone
            else:
                logger.error(f"A JKTS monitor (pid {old_pid}) still appears to be running in this directory. "
                             f"Stop it first (kill {old_pid}) or clear 'monitor_pid' in .metadata if it is stale.")
                sys.exit(1)
    update_metadata(runtime.start_dir, monitor_pid=os.getpid())


def restore_settings(parser):
    meta = load_metadata(runtime.start_dir)
    if not meta:
        return
    restored = []

    # Reaction type: -OH/-Cl/-NO3 are store_true flags; none given means all False
    reaction = meta.get('reaction')
    if reaction in ('OH', 'Cl', 'NO3') and not (runtime.args.OH or runtime.args.Cl or runtime.args.NO3):
        setattr(runtime.args, reaction, True)
        restored.append(f"reaction={reaction}")

    # QC backend: -G16/-ORCA are store_true flags; none given means all False
    program = meta.get('program')
    if program in ('G16', 'ORCA') and not (runtime.args.G16 or runtime.args.ORCA):
        runtime.args.ORCA = program == 'ORCA'
        runtime.args.G16 = program == 'G16'
        runtime.QC_program = program
        restored.append(f"program={program}")

    slurm = meta.get('slurm', {})
    for key, value in (('method', meta.get('method')),
                       ('basis_set', meta.get('basis_set')),
                       ('T', meta.get('temperature')),
                       ('par', slurm.get('par')),
                       ('time', slurm.get('time')),
                       ('cpu', slurm.get('cpu')),
                       ('mem', slurm.get('mem'))):
        if value is not None and getattr(runtime.args, key) == parser.get_default(key):
            setattr(runtime.args, key, value)
            restored.append(f"{key}={value}")

    if restored:
        logger.info(f"Restored settings from .metadata: {', '.join(restored)}")


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
