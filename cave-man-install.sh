#!/bin/sh -e

##########################################################################
#
#   Perform a cave-man install for development and testing purposes.
#   For production use, this software should be installed via a package
#   manager such as Debian packages, FreeBSD ports, MacPorts, pkgsrc, etc.
#       
#   History:
#   Date        Name        Modification
#   2021-07-12  Jason Bacon Begin
##########################################################################

# Default to ../local if PREFIX is not set
: ${PREFIX:=../local}

# OS-dependent tricks
# Set rpath to avoid picking up libs installed by package managers in
# /usr/local/lib, etc.
case $(uname) in
*)
    if [ -z "$CFLAGS" ]; then
	export CFLAGS="-Wall -g -O"
    fi
    LIBDIR=$(realpath $PREFIX/lib)
    export LDFLAGS="-L$LIBDIR -Wl,-rpath,$LIBDIR:/usr/lib:/lib"
    for pkgsrc in /usr/pkg /opt/pkg ~/Pkgsrc/pkg; do
	if [ -e $pkgsrc ]; then
	    echo "Using $pkgsrc..."
	    export LOCALBASE=$pkgsrc
	fi
    done
    ;;

esac

make clean
make depend
make install
