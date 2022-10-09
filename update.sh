#!/usr/bin/env bash

set -e
set -x

SRC="tmp"
DEST="src/parmetis"

git clone https://github.com/KarypisLab/ParMETIS.git $SRC
mkdir -p $DEST
mkdir -p $DEST/libparmetis
mkdir -p $DEST/include

cp -f $SRC/LICENSE $DEST/LICENSE
cp -fa $SRC/include/*.h $DEST/include/
cp -fa $SRC/libparmetis/*.{h,c} $DEST/libparmetis/
# cp -fa $SRC/programs/*.{h,c} app/
cp -fa $SRC/graphs ./

# Remove files that are not needed
rm -rf $SRC
