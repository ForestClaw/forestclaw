#! /bin/bash

# Build and install the latest ForestClaw develop branch including p4est,
# libsc, zlib and jansson.  We first download and install a recent zlib to
# a local directory, then do the same for the jansson library.  Then we
# clone and install the current develop branches of p4est and libsc using
# that zlib and jansson installation and finally clone and install the
# current develop branch of ForestClaw, linking against all of the above.

# This results in five installation directories that any higher
# level software package may be compiled and linked against.
# The options are similar to those used in this script.
# In particular, the -rpath option may turn out useful.

# let the p4est installation script set the stage
BLP4="p4est-build-wdeps.sh"
wget -N "https://github.com/cburstedde/p4est/raw/develop/doc/$BLP4"     && \
source "$BLP4"                                          && \
rm "$BLP4"                                              || exit 1

# provide environment for the Fortran compile
export LIBS="$LIBS -lgfortran"

# clone, build and install ForestClaw
git clone --depth 1 https://github.com/forestclaw/forestclaw.git -b develop && \
cd forestclaw                                           && \
./bootstrap "$PREFIX/libsc/share/aclocal" "$PREFIX/p4est/share/aclocal" && \
mkdir build                                             && \
cd build                                                && \
../configure $CONFIG --with-sc="$PREFIX/libsc" --with-p4est="$PREFIX/p4est" \
                     --prefix="$PREFIX/fclaw"           && \
(make -j V=0 || make -j V=0 || make -j V=0)             && \
make -j install V=0                                     && \
cd ../../                                               && \
rm -rf forestclaw/.git                                  && \
rm -r forestclaw                                        || bdie "ForestClaw"
