#!/bin/bash
source /usr/local/gromacs/bin/GMXRC
# This is designed to restart an equilibration, either NVT or NPT, run.
# Script to run a single gmx mdrun job and restart it if it fails
# This script will move right into a restart without any user input
# if a failure is detected.
# This will attempt a maximum number of restarts as specified in
# the while loop control.

# Gives a more informative error when something goes wrong
# with the script.
error_report() {
echo "Error $1 on j $2, iMCMC $3, stage $4"
exit 1
}

# Assemble the parameters required for the run and potential restart
clean() {   # Make it so that everything is killed on an interrupt
local pids=$(jobs -pr)
echo "On exit sending kill signal to: $pids"
[ -n "$pids" ] && kill $pids
exit 1
}
trap "clean" SIGINT SIGTERM EXIT SIGQUIT  # Call clean for any such signal

cass_path="$1"
inp_file_2NVT="$2"
inp_file_equil="$3"
inp_file_prod="$4"

"$cass_path" "$inp_file_2NVT" > run_info_2NVT 2>&1  # Run 2NVT Cassandra
"$cass_path" "$inp_file_equil" > run_info_equil 2>&1  # Run GEMC equil Cassandra
"$cass_path" "$inp_file_prod" > run_info_prod 2>&1 # Run GEMC prod Cassandra
