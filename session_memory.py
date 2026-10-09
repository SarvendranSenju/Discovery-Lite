"""
Small JSON-backed sidecar that remembers recently used receptor/ligand paths
and PDB IDs between GUI sessions, so fields like "8FK4" or a reused
native_ligand.sdf path don't have to be retyped every time the app reopens.

The file lives next to the GUI scripts as session_memory.json:
    {
      "last_used": {"redock_pdb": "8FK4", "dock_receptor": "...", ...},
      "recent":    {"redock_pdb": ["8FK4", "3AIC", ...], ...}
    }

It's plain JSON, so it can be inspected or hand-edited, or just deleted to
reset. Writes are atomic (write to a .tmp file then os.replace) so a crash
mid-write can't corrupt it.
"""

import json
import os
import threading

MAX_RECENT = 12


class SessionMemory:
    def __init__(self, path):
        self.path = path
        self._lock = threading.Lock()
        self._data = {"last_used": {}, "recent": {}}
        self.load()

    def load(self):
        if not os.path.isfile(self.path):
            return
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
            if isinstance(loaded, dict):
                self._data["last_used"] = loaded.get("last_used") or {}
                self._data["recent"] = loaded.get("recent") or {}
        except (json.JSONDecodeError, OSError):
            # Corrupt or unreadable sidecar - start fresh rather than crash the GUI.
            pass

    def save(self):
        with self._lock:
            try:
                tmp_path = self.path + ".tmp"
                with open(tmp_path, "w", encoding="utf-8") as f:
                    json.dump(self._data, f, indent=2, sort_keys=True)
                os.replace(tmp_path, self.path)
            except OSError:
                pass

    def remember(self, key, value):
        """Record `value` as the last-used and most-recent entry for `key`."""
        value = (value or "").strip()
        if not value:
            return
        self._data["last_used"][key] = value

        recent = self._data["recent"].setdefault(key, [])
        recent[:] = [v for v in recent if v.lower() != value.lower()]
        recent.insert(0, value)
        del recent[MAX_RECENT:]

        self.save()

    def get_last(self, key, default=""):
        return self._data["last_used"].get(key, default)

    def get_recent(self, key):
        return list(self._data["recent"].get(key, []))
