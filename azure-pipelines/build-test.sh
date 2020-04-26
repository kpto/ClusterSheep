#!/bin/bash
# Based on https://github.com/scikit-hep/azure-wheel-helpers/blob/master/build-wheels.sh

# Collect the pythons
pys=(/opt/python/*/bin)

# Filter out Python 2
pys=(${pys[@]//*2*/})

# Print list of Python's being used
echo "Using Pythons: ${pys[@]}"

# Compile wheels
for PYBIN in "${pys[@]}"; do
    "${PYBIN}/pip" install -r /io/$build_requirements_file
    "${PYBIN}/pip" wheel -v /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/$package_name-*.whl; do
    auditwheel repair --plat $PLAT "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in "${pys[@]}"; do
    VERSION=$($PYBIN/python -V)
    "${PYBIN}/pip" install -r /io/$test_requirements_file
    "${PYBIN}/pip" install $package_name --no-index -f /io/wheelhouse
    "${PYBIN}/pytest" --nunitxml="/io/${VERSION// /}-test-output.xml" /io/tests
done