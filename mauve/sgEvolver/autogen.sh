#!/bin/sh
mkdir -p config
autoheader
autoreconf --force --install -I config && \
echo "Now run ./configure --prefix=$HOME ; make install"

