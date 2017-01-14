#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <platformtag>"
    exit 2
fi

platformtag="$1"

set -eu

ippdir=/opt/intel/ipp

echo
if [ -d "$ippdir" ]; then
    echo "Found IPP directory $ippdir, considering IPP as well as other options"
else
    echo "No IPP directory $ippdir, not testing with IPP"
fi

tmpfile=$(mktemp -p "" "test_XXXXXX")

run() {
    successtext="$1"
    shift
    echo -n "Running \"$@\"..."
    if "$@" > "$tmpfile" 2>&1 ; then
	if fgrep -q "$successtext" "$tmpfile" ; then
	    echo " OK"
	    return 0
	else
	    echo " Failed"
	    cat "$tmpfile"
	    return 1
	fi
    else
	echo " Failed"
	cat "$tmpfile"
	return 1
    fi
}

for mf in Makefile build/Makefile.$platformtag build/Makefile.$platformtag.* ; do

    case "$mf" in
	*ipp)
	    if [ ! -d "$ippdir" ]; then
		continue
	    fi;;
    esac

    echo
    echo "Building and testing with $mf:"
    echo
    
    make -f "$mf" clean >/dev/null
    run "No errors detected" make -f "$mf" test

    for t in test-* ; do
	if [ -x "$t" ]; then
	    run "no leaks are possible" valgrind --leak-check=full ./"$t"
	fi
    done
done


