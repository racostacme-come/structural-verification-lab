import json
import subprocess
import sys


def test_cli_json_and_error(tmp_path):
    path = tmp_path / "nested" / "result.json"
    command = [sys.executable, "-m", "bucklemap.cli", "solve"]
    subprocess.run([*command, "--spring", "10", "--output", str(path)], check=True)
    data = json.loads(path.read_text())
    assert 28 < data["loads"][0] < 29
    bad = subprocess.run([*command, "--elements", "0"], capture_output=True, text=True)
    assert bad.returncode == 2 and "elements must be >= 1" in bad.stderr
