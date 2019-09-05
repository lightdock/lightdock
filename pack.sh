#!/bin/bash

# Remove old builds
rm -rf build dist lightdock.egg-info/

# Build dist
python3 setup.py sdist bdist_wheel

# Deploy on the testing server
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

