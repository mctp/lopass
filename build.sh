#!/usr/bin/env bash
LOPASS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "working in: "$LOPASS_DIR
mkdir -p $LOPASS_DIR/boost

## build htslib
cd $LOPASS_DIR
wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar -xf htslib-1.12.tar.bz2
rm htslib-1.12.tar.bz2
mv htslib-1.12 htslib
cd htslib
make -j4

## build boost
cd $LOPASS_DIR/boost
wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2
tar --bzip2 -xf boost_1_73_0.tar.bz2
rm boost_1_73_0.tar.bz2
mv boost_1_73_0 src
cd src
./bootstrap.sh --with-libraries=iostreams,program_options --prefix=../build
./b2 install

## patching glimpse
cd $LOPASS_DIR/GLIMPSE
BOOST_LIB=$LOPASS_DIR/boost/build/lib/
for i in */makefile; do
    echo "patching: " $i
    sed -i.bkp 's,'"/usr/local/lib/"','"$BOOST_LIB"',g' $i
done

## building glimpse
cd $LOPASS_DIR/GLIMPSE
make -j4

