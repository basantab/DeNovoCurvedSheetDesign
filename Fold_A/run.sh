#!/bin/bash
# Argument 1 is your rosetta_scripts executable path
# Argument 2 is the max runtime in seconds
for i in $(ls -d r*-all_rev*/); do
  cd $i;
  mkdir run_dir;
  cd run_dir;
  $1 -s ../input.pdb -parser:protocol ../input.xml -nstruct 10000 -maxruntime $2 @../flags;
  cd ..;
done
