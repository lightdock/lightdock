#!/bin/bash

## Compiles all the required Python extensions found in the LightDock project

# Displays a coloured message
function cinfo() {
    COLOR='\033[01;33m' # bold yellow
    RESET='\033[00;00m' # normal white
    MESSAGE=${@:-"${RESET}No message passed"}
    echo -e "${COLOR}${MESSAGE}${RESET}"
}

cinfo "Compiling extensions..."
CURRENT_PATH=`pwd`

for i in `find . -name "compile.sh"`;do
	dest=`dirname ${i}`;
	cinfo $dest;
	cd $dest;
	./compile.sh;
	cd $CURRENT_PATH;
done

cinfo "Done."
