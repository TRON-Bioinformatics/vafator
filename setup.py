from setuptools import find_packages, setup
import vafator


VERSION = vafator.VERSION

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup()