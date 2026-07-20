import sys

from pufscan.cli import app

if __name__ == "__main__":
    sys.argv.insert(1, "prepare-gtex")
    app()

