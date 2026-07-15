import os
import pickle
from threading import Lock

from output import logger
import runtime

# Serializes checkpoint writes across the monitoring threads of one process.
_checkpoint_lock = Lock()


def checkpoint_path(directory):
    base = os.path.basename(os.path.normpath(directory))
    return os.path.join(directory, f"{base}_checkpoint.pkl")


def atomic_pickle_dump(obj, file_path):
    tmp_path = file_path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump(obj, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, file_path)


def load_pickle(file_path):
    try:
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    except (FileNotFoundError, pickle.UnpicklingError, EOFError, AttributeError, OSError):
        return None


def save_checkpoint(molecules):
    if not molecules:
        return
    directory = molecules[0].directory or runtime.start_dir
    if not directory:
        return
    with _checkpoint_lock:
        with runtime.global_molecules_lock:
            payload = list(molecules) + [g for g in runtime.global_molecules if g not in molecules]
        try:
            atomic_pickle_dump(payload, checkpoint_path(directory))
        except (OSError, pickle.PicklingError) as e:
            logger.warning(f"Could not write checkpoint in {directory}: {e}")


def load_checkpoint(directory):
    return load_pickle(checkpoint_path(directory))
