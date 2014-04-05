#!/bin/bash
echo "run from root of repository checkout!"
git fetch ngspice master
git subtree pull --prefix=s/ext/ngspice ngspice master --squash

