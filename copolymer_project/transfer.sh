#!/bin/bash

cd replicates
for dir in */; do
	if [ -d "$dir" ]; then

		cd "$dir"
		cp *.png ../../transfer
		cd ..
	fi
done
cd ..
