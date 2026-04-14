import csv
import subprocess
import sys

MANIFEST = "/opt/tools/build_manifest.csv"

def get_version(tool):
    try:
        out = subprocess.check_output([tool, "--version"], stderr=subprocess.STDOUT)
        return out.decode().strip()
    except Exception:
        return None

def main():
    with open(MANIFEST) as f:
        reader = csv.DictReader(f)

        for row in reader:
            tool = row["tool"]
            expected_version = row["version"]

            actual = get_version(tool)

            if actual is None:
                print(f"[FAIL] {tool} not found")
                sys.exit(1)

            if expected_version not in actual:
                print(f"[FAIL] {tool}")
                print(f"expected: {expected_version}")
                print(f"got:      {actual}")
                sys.exit(1)

            print(f"[OK] {tool} {actual}")

if __name__ == "__main__":
    main()