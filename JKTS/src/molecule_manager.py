"""Thread-safe molecule manager used as a drop-in replacement for a global list.

Expose a module-level instance `global_molecules` so existing code that
references that name can keep working with minimal changes.
"""
from threading import RLock
import pickle
from typing import Iterable, List, Iterator, Optional, Callable


class MoleculeManager:
    def __init__(self, molecules: Iterable = ()):
        self._lock = RLock()
        self._molecules: List = list(molecules)

    # List-like mutators
    def append(self, mol) -> None:
        with self._lock:
            self._molecules.append(mol)

    def extend(self, mols: Iterable) -> None:
        with self._lock:
            self._molecules.extend(list(mols))

    def remove(self, mol) -> None:
        with self._lock:
            self._molecules.remove(mol)

    def clear(self) -> None:
        with self._lock:
            self._molecules.clear()

    # List-like accessors
    def as_list(self) -> List:
        with self._lock:
            return list(self._molecules)

    def to_list(self) -> List:
        return self.as_list()

    def __len__(self) -> int:
        with self._lock:
            return len(self._molecules)

    def __bool__(self) -> bool:
        return len(self) > 0

    def __iter__(self) -> Iterator:
        # Return an iterator over a snapshot to avoid concurrent modification issues
        return iter(self.as_list())

    def __getitem__(self, idx):
        with self._lock:
            return self._molecules[idx]

    # Convenience helpers
    def find(self, predicate: Callable) -> Optional:
        with self._lock:
            for m in self._molecules:
                if predicate(m):
                    return m
            return None

    def filter(self, predicate: Callable) -> List:
        with self._lock:
            return [m for m in self._molecules if predicate(m)]

    def group_by(self, keyfunc: Callable) -> dict:
        with self._lock:
            groups = {}
            for m in self._molecules:
                k = keyfunc(m)
                groups.setdefault(k, []).append(m)
            return groups

    def to_pickle(self, path: str) -> None:
        with self._lock:
            with open(path, 'wb') as f:
                pickle.dump(self._molecules, f)


# Module-level instance to be used by other modules that expect a
# `global_molecules` list. This keeps the migration minimal.
global_molecules = MoleculeManager()
