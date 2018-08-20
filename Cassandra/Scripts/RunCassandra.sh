#!/bin/bash
# Run Cassandra simulations in sets.

# Gives a more informative error when something goes wrong
# with the script.
#
##########################################################
# USE INSTRUCTIONS:
# Launch this script as ./RunCassandra.sh & then run top
# to prevent the session from timing out.
#

#########################################################
# HELPER FUNCTIONS
error_report() {
echo "Error in $PWD: $1"
exit 1
}

clean() {   # Make it so everything is killed on an interrupt
local pids=$(jobs -pr)
echo "On exit from $PWD sending kill signal to: $pids"
[ -n "$pids" ] && kill $pids
exit 1
}
trap "clean" SIGINT SIGTERM EXIT SIGQUIT  # Call cleanup when asked to
##############################################################

######### INPUT SPECIFICATIONS ##############################
scripts_dir="$PWD"
working_dir="$PWD"/Sequential_trial3_highTemp
cass_path="/home/ram9/Cassandra_Helium/Src/cassandra_gfortran.exe"
inp_file_2NVT="helium.2NVT.inp" # File to actually run
inp_file_equil="helium.equil.inp" # File to actually run
inp_file_prod="helium.prod.inp" # File to actually run
inp_dir="$PWD"/input  # Contains all files necessary for a run
# Each temperature set requires a step number specification
steps_NVT=(100000 100000 100000 100000 100000)
steps_equil=(1000000 1000000 1000000 1000000 1000000)
steps_prod=(5000000 5000000 5000000 5000000 5000000)
#temps1=(7 8 9 10 11)  # Temperatures in each box
#temps2=(7 8 9 10 11)
# What different size conditions for each temperature
#sizes1=(2400 2400 2400 2400 2400)  # Molecules in each box for each temperature
#sizes2=(400 400 400 400 400)
#boxes1=(37.40893 38.04369 38.78536 39.70215 40.99157)  # Box sizes for each temperature
#boxes2=(155.03813 107.04542 82.79921 62.73321 51.69863)
temps1=(11.5 12 12.5 13 13.5)  # Temperatures in each box
temps2=(11.5 12 12.5 13 13.5)
# What different size conditions for each temperature
sizes1=(2400 2400 2400 2400 2400)  # Molecules in each box for each temperature
sizes2=(400 400 400 400 400)
boxes1=(40.99157 40.99157 40.99157 40.99157 40.99157)  # Box sizes for each temperature
boxes2=(51.69863 51.69863 51.69863 51.69863 51.69863)

# What different size conditions for each temperature
#sizes1=(1200 1200 1200 1200 1200)  # Molecules in each box for each temperature
#sizes2=(200 200 200 200 200)
#boxes1=(29.69148742 30.19529676 30.78396063 31.51161734 32.53503067)  # Box sizes for each temperature
#boxes2=(123.0538453 84.96200616 65.71777653 49.79138177 41.03322982)
# How many replicates of each temperature and size combination
replicates=1
# How many processes should be allowed to run at any time
max_jobs=5
# Random seeds to start with; they will be incremented for each job
seed1=1211131639
seed2=1211131640

#################################################################
################### DIRECTORY PREPARATION ######################


cd "$working_dir" || error_report "Cannot change to $working_dir"
# Print start up informaiton
echo "Running ${#temps1[@]} temperatures at ${#sizes1[@]}"
echo "states each with $replicates replicates per state"
echo "Temperatures 1: ${temps1[@]}"
echo "Temperatures 2: ${temps2[@]}"
echo "Step number at each temperature: ${steps[@]}"
echo "Molecule numbers 1: ${sizes1[@]}"
echo "Molecule numbers 2: ${sizes2[@]}"
echo "Boxes 1: ${boxes1[@]}"
echo "Boxes 2: ${boxes2[@]}"

# Adjust for ` indexing
temp_nums=${#temps1[@]}
temp_nums=$((temp_nums - 1))
rep_nums=$((replicates - 1))
size_nums=${#sizes1[@]}
size_nums=$((size_nums - 1))


for i in $(seq 0 $temp_nums); do
  mkdir "temp_$i"
  cd "temp_$i" || error_report "Cannot change to temp_$i"
  for j in $(seq 0 $rep_nums); do
    cp -r "$inp_dir" "rep_$j"  # Copy input files
    cd "rep_$j" || error_report "Cannot change to rep_$j"
    # Copy in the input file and the fragment directory
    #cp "$working_dir"/"$inp_file" "$inp_file"
    #cp "$working_dir"/"$frag_dir" "$frag_dir"    
    # Subsitute in temperatures, steps, replicates, sizes
 
    for k in $(seq 0 2); do

      if [ $k -eq 0 ]; then
         
        inp_file="$inp_file_2NVT"
        steps=("${steps_NVT[@]}")  # Need [@] to get entire array, need {} as usual, and need () to create new array

      elif [ $k -eq 1 ]; then

        inp_file="$inp_file_equil"
        steps=("${steps_equil[@]}")

      elif [ $k -eq 2 ]; then

        inp_file="$inp_file_prod"
        steps=("${steps_prod[@]}")

      fi

      sed -i -e s/insert_box_1/${boxes1[i]}/ "$inp_file"
      sed -i -e s/insert_box_2/${boxes2[i]}/ "$inp_file"
      sed -i -e s/insert_temps_1/${temps1[i]}/ "$inp_file"
      sed -i -e s/insert_temps_2/${temps2[i]}/ "$inp_file"
      sed -i -e s/insert_mols_1/${sizes1[i]}/ "$inp_file"
      sed -i -e s/insert_mols_2/${sizes2[i]}/ "$inp_file"
      sed -i -e s/insert_steps/${steps[i]}/ "$inp_file"
      sed -i -e s/insert_seed_1/$seed1/ "$inp_file"
      sed -i -e s/insert_seed_2/$seed2/ "$inp_file"

    done
    
    seed1=$((seed1 + 1))
    seed2=$((seed2 + 1))
    cd .. || error_report "Cannot change to .."
  done
  cd .. || error_report "Cannot change to .."
done

######################################################

################ CASSANDRA RUNNING ####################

cur_jobs=0
# For each temperature
for i in $(seq 0 $temp_nums); do
  cd "temp_$i" || error_report "Cannot change to temp_$i"
  # For each replicate
  for j in $(seq 0 $rep_nums); do
    cd "rep_$j" || error_report "Cannot change to rep_$j"
    #"$cass_path" "$inp_file" > run_info 2>&1 &  # Start Cassandra
    bash "$scripts_dir"/Run2NVT_GEMCequil_GEMCprod.sh "$cass_path" "$inp_file_2NVT" "$inp_file_equil" "$inp_file_prod"  & #Run the sequential 2NVT, GEMC equil, and GEMC prod
    cur_jobs=$((cur_jobs + 1))
    # echo "cur_jobs $cur_jobs"
    #if [ $cur_jobs -ge $max_jobs ]; then
      # Wait for the jobs to finish.
      # echo "Waiting for jobs batch (cur_jobs $cur_jobs)"
      # Wait for all those processes to finish
    #  wait || error_report "Failed to wait for all jobs!"
      # echo "Waited: processes running: "
      # ps -a
    #  cur_jobs=0  # Start the next set
    #fi
    cd .. || error_report "Cannot change to .."
  done
  cd .. || error_report "Cannot change to .."
done

top

# Wait for any straggler jobs
wait

#########################################################
### I Felt the need to make this symmetrical ############
