#!/bin/sh

# set -x

if [ -z "$SPICE_SCRIPTS" ] ; then
    SPICE_SCRIPTS=`dirname $0`
    export SPICE_SCRIPTS
    if [ -z "$ngspice_vpath" ] ; then
        ngspice_vpath=.
        export ngspice_vpath
    fi
fi

# ls -ld $(realpath $SPICE_SCRIPTS) $SPICE_SCRIPTS/spinit
# echo "---ngspice_vpath = $ngspice_vpath"

SPICE=$1
TEST=$2

FILTER="CPU|Dynamic|Note|Circuit|Trying|Reference|Date|Doing|---|v-sweep|time|est|Error|Warning|Data|Index|trans|acan|oise|nalysis|ole|Total|memory|urrent|Got|Added|BSIM|bsim|B4SOI|b4soi|codemodel|^binary raw file|^ngspice.*done"

testname=`basename $TEST .cir`
testdir=`dirname $TEST`

HOST_TYPE=`uname -srvm`

case $HOST_TYPE in
    MINGW32*)
      $SPICE --batch $testdir/$testname.cir -o $testname.test &&\
      sed -e 's/e-000/e+000/g' $testname.test | sed 's/e-0/e-/g' | sed 's/e+0/e+/g' | egrep -v "$FILTER" > $testname.test_tmp &&\
      sed -e 's/-0$/ 0/g' $testdir/$testname.out | egrep -v "$FILTER" > $testname.out_tmp
      if diff -B -w -u $testname.out_tmp $testname.test_tmp; then
          rm $testname.test $testname.test_tmp $testname.out_tmp
          exit 0
      fi
      rm -f $testname.test_tmp $testname.out_tmp
      sed -e 's/e-000/e+000/g' $testname.test | sed 's/e-0/e-/g' | sed 's/e+0/e+/g' > $testname.test_tmp
      mv $testname.test_tmp $testname.test
      ;;
    Linux*|Darwin*|CYGWIN*)
      $SPICE --batch $testdir/$testname.cir >$testname.test &&\
      egrep -v "$FILTER" $testname.test > $testname.test_tmp &&\
      egrep -v "$FILTER" $testdir/$testname.out > $testname.out_tmp
      if diff -B -w -u $testname.out_tmp $testname.test_tmp; then
          rm $testname.test $testname.test_tmp $testname.out_tmp
          exit 0
      fi
      rm -f $testname.test_tmp $testname.out_tmp
      ;;
    SunOS*|OpenBSD*)
      $SPICE --batch $testdir/$testname.cir >$testname.test &&\
      sed -e '/^$/d' $testname.test | egrep -v "$FILTER" > $testname.test_tmp &&\
      sed -e '/^$/d' $testdir/$testname.out | egrep -v "$FILTER" > $testname.out_tmp
      if diff -b -w $testname.out_tmp $testname.test_tmp; then
          rm $testname.test $testname.test_tmp $testname.out_tmp
          exit 0
      fi
      rm -f $testname.test_tmp $testname.out_tmp
      ;;
    *)
      echo Unknown system type!
      echo $HOST_TYPE
      echo ./tests/bin/checks.sh may need updating for your system
      ;;
esac

exit 1
