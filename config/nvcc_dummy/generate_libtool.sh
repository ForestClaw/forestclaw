#!/bin/sh
# first argument is path to nvcc, second argument is path to root of source
echo "Putting nvcclibtool in ${PWD} from files in ${2}"
build_dir=$PWD
cd $2/config/nvcc_dummy
aclocal
autoconf
CC=$1 ./configure
cp config.log ${build_dir}/nvcc_config.log
cp libtool ${build_dir}/nvcclibtool
rm -rf configure libtool config.log config.status autom4te.cache aclocal.m4
