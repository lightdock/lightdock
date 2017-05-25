#!/bin/bash

# Displays a coloured message
function cinfo() {
    COLOR='\033[01;33m' # bold yellow
    RESET='\033[00;00m' # normal white
    MESSAGE=${@:-"${RESET}No message passed"}
    echo -e "${COLOR}${MESSAGE}${RESET}"
}

LIB_TEST=true
REG_TEST=true
if [ "$#" -ne 1 ]
then
  echo "Usage: $0 [lib|reg|all]"
  exit 1
fi

if [ "$1" = "lib" ] ; then
    REG_TEST=false
fi

if [ "$1" = "reg" ] ; then
    LIB_TEST=false
fi

if [ "$LIB_TEST" = true ] ; then
    cinfo "******* Library tests *******"
    cd lightdock
    nosetests-2.7
    echo ""
    cd ..
fi

if [ "$REG_TEST" = true ] ; then
    cinfo "******* Regression tests *******"
    cd bin
    nosetests-2.7
    cd ..
    cinfo "Done."
fi
