#!/bin/bash -e

# Set up an environment to run tests under Travis CI (see ../.travis.yml)

if [ $# -ne 1 ]; then
  echo "Usage: $0 python_version"
  exit 1
fi

cur_dir=$(pwd)
python_version=$1
temp_dir=$(mktemp -d)

cd ${temp_dir}

if [ ${python_version} = "2.7" ]; then
  BOOST=""
  pip="pip<=19.3.1"
  # Python.h includes crypt.h, which is no longer provided by default
  crypt="libxcrypt"
else
  BOOST="libboost-devel"
  pip="pip"
  crypt=""
fi
conda config --remove channels defaults  # get conda-forge, not main, packages
conda create --yes -q -n python${python_version} -c salilab -c conda-forge python=${python_version} ${pip} ${crypt} scipy matplotlib imp-nightly ${BOOST} gxx_linux-64 eigen cereal swig cmake
eval "$(conda shell.bash hook)"
conda activate python${python_version}
pip install pytest-cov coverage

cd ${cur_dir}

rm -rf ${temp_dir}
