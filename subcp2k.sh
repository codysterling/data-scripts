#!/bin/bash

# Colors
normal=`tput sgr0`
bold=`tput bold`
red=`tput setaf 1`
yellow=`tput setaf 3`
green=`tput setaf 2`
cyan=`tput setaf 6`
magenta=`tput setaf 5`
blue=`tput setaf 4`
underline=`tput smul`

# Explaining code:
if [[ $1 == "" ]]
then
  echo "Usage:"
  echo "subcp2k job(.inp) -option arg ..."
  echo "${underline}Options:${normal}"
  echo "${magenta}-N: nodes: integer number of nodes to use (default 2)${normal}"
  echo "${bold}${green}-t: walltime: #d (days 1-15), # (hrs 1-360), or days-hrs:min:sec (default 1-00:00:00)${normal}"
  exit
fi

# Stripping .chm from input if there
job=$(echo ${1%.*})

# Checking that input file is there
if [[ ! -f $job.inp ]]
then
  echo "$1 does not exist!"
  exit
fi

# Putting subchem arguments into an array so arguments can be pulled generally
for i in $@
do
  args+=($i)
done

# This function finds the indices of the argument in the array args
# Sourced from http://stackoverflow.com/questions/30770780/find-the-position-of-an-element-in-a-list-in-bash
indexof() { i=-1; for ((j=0;j<${#args[@]};j++)); do [ "${args[$j]}" = "$1" ] && { i=$j; break; } done; echo $i; }

# Sourcing input file and adding job name
jobscript=cp2k-submit.sh
cp /home/x_codst/proj/cp2k-templates/$jobscript .
sed -i "s:JOBNAME:$job:g" ./$jobscript

# Checking other flags, setting options as needed
optflag=`indexof "-q"`
if [[ $optflag != "-1" ]]
then
  argindex=$(($optflag+1))
  argval=${args[$argindex]}
  sed -i "s:-p normal:-p $argval:" ./$jobscript
  echo "${bold}${blue}Submitted to $argval queue${normal}"
fi

optflag=`indexof "-N"`
if [[ $optflag != "-1" ]]
then
  argindex=$(($optflag+1))
  argval=${args[$argindex]}
  sed -i "s:-N 2:-N $argval:" ./$jobscript
  sed -i "s:-n 64:-n $(($argval*32)):" ./$jobscript
  echo "${magenta}Submitted with $argval nodes${normal}"
fi

optflag=`indexof "-c"`
if [[ $optflag != "-1" ]]
then
  argindex=$(($optflag+1))
  argval=${args[$argindex]}
  sed -i "s:-n 32:-n $argval:" ./$jobscript
  echo "${bold}${cyan}Submitted with $argval cores${normal}"
fi

optflag=`indexof "-t"`
if [[ $optflag != "-1" ]]
then
  argindex=$(($optflag+1))
  argval=${args[$argindex]}
  if [[ $argval =~ ^-?[0-9]+d+$ ]]
  then
    argval=$(echo ${argval%d})
    if [[ $argval -le 30 ]]
    then
      sed -i "s/-t 1-00:00:00/-t $argval-00:00:00/g" ./$jobscript
      echo "${bold}${green}Submitted with $argval days${normal}"
    else
      echo "Your time of $argval days is out of range, not submitted."
      exit
    fi
  elif [[ $argval =~ ^-?[0-9]+$ ]]
  then
    if [[ $argval -le 360 ]]
    then
      sed -i "s/-t 1-00:00:00/-t $argval:00:00/g" ./$jobscript
      echo "${bold}${green}Submitted with $argval hours${normal}"
    else
      echo "Your time of $argval hours is out of range, not submitted."
      exit
    fi
  else
    sed -i "s/-t 1-00:00:00/-t $argval/g" ./$jobscript
    echo "${bold}${green}Submitted with time $argval${normal}"
  fi
fi

# Submitting job
# Here using the qsub SLURM wrapper to allow job.eX and job.oX files
echo "${bold}${magenta}Job ID: `sbatch -J $job $jobscript`${normal}"
