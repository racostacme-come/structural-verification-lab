"""The provenance gate must reject changed or missing snapshot bytes."""
import hashlib

import pytest

import suite


def test_snapshot_integrity_gate(tmp_path, monkeypatch):
    path = tmp_path / "modules" / "sample" / "source.cpp"
    path.parent.mkdir(parents=True)
    original = b"int main() { return 0; }\n"
    path.write_bytes(original)
    monkeypatch.setattr(suite, "ROOT", tmp_path)
    monkeypatch.setattr(suite, "MANIFEST", {"modules": [{
        "id": "sample", "sha256": {"source.cpp": hashlib.sha256(original).hexdigest()}
    }]})
    suite.verify()
    path.write_bytes(original.replace(b"return 0", b"return 1"))
    with pytest.raises(RuntimeError, match="Imported source changed"):
        suite.verify()
    path.unlink()
    with pytest.raises(FileNotFoundError):
        suite.verify()
