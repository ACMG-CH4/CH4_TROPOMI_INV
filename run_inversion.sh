#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH -t 0-03:00
#SBATCH --mem 4000
#SBATCH -o run_inversion_%j.out
#SBATCH -e run_inversion_%j.err

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Function to parse yaml files from shell script
# By Stefan Farestam via stackoverflow:
# https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Get configuration
eval $(parse_yaml ../../../setup_CH4/config.yml)
# This defines $PriorError, $ObsError, and $Gamma (in addition to the other config entries)
# Parsing the config file here facilitates generation of inversion ensembles
# All that needs to be done is to edit the config file for $PriorError, $ObsError, and $Gamma, then re-run this script

#=======================================================================
# Configuration (these settings generated on initial setup)
#=======================================================================
StartDate={START}
EndDate={END}
LonMin={LON_MIN}
LonMax={LON_MAX}
LatMin={LAT_MIN}
LatMax={LAT_MAX}
nElements={STATE_VECTOR_ELEMENTS}
nBufferClusters={BUFFER_CLUSTERS}
MyPath={MY_PATH}
RunName={RUN_NAME}
isAWS={IS_AWS}
SpinupDir="${MyPath}/${RunName}/spinup_run"
JacobianRunsDir="${MyPath}/${RunName}/jacobian_runs"
PosteriorRunDir="${MyPath}/${RunName}/posterior_run"
StateVectorFile={STATE_VECTOR_PATH}
SensiDir="./Sensi"
GCDir="./data_GC"
JacobianDir="./data_converted"
TROPOMIDir="./data_TROPOMI"

# Download TROPOMI data. Disable this setting if rerunning for a time period
# to avoid redownloading existing data
if "$isAWS"; then
    FetchTROPOMI=true
else
    FetchTROPOMI=false
fi
    
# Only matters for Kalman filter multi-week inversions
FirstSimSwitch=true

# Load the Python environment
# Comment this out for now and instruct users to activate their own Python envs
# module load Anaconda3/5.0.1-fasrc01
# source activate tropomi_oilgas_jupyter_3.6


printf "\n=== EXECUTING RUN_INVERSION.SH ===\n"
    
#=======================================================================
# Error checks
#=======================================================================

# Make sure specified paths exist
if [[ ! -d ${JacobianRunsDir} ]]; then
    printf "${JacobianRunsDir} does not exist. Please fix JacobianRunsDir in run_inversion.sh.\n"
    exit 1
fi
if [[ ! -f ${StateVectorFile} ]]; then
    printf "${StateVectorFile} does not exist. Please fix StateVectorFile in run_inversion.sh.\n"
    exit 1
fi

#=======================================================================
# Postprocess the SpeciesConc and LevelEdgeDiags files from GEOS-Chem
#=======================================================================

printf "Calling postproc_diags.py, FSS=$FirstSimSwitch\n"
if "$FirstSimSwitch"; then
    if [[ ! -d ${SpinupDir} ]]; then
	printf "${SpinupDir} does not exist. Please fix SpinupDir or set FirstSimSwitch to False in run_inversion.sh.\n"
	exit 1
    fi
    PrevDir=$SpinupDir
else
    PrevDir=%$PosteriorRunDir
    if [[ ! -d ${PosteriorRunDir} ]]; then
	printf "${PosteriorRunDir} does not exist. Please fix PosteriorRunDir in run_inversion.sh.\n"
	exit 1
    fi
fi
printf "  - Hour 0 for ${StartDate} will be obtained from ${PrevDir}\n"

python postproc_diags.py $RunName $JacobianRunsDir $PrevDir $StartDate; wait
printf "DONE -- postproc_diags.py\n\n"

#=======================================================================
# Calculate GEOS-Chem sensitivities and save to Sensi directory
#=======================================================================

# 50% perturbation implied by PerturbValue (config.yml): 1.5 - 1 = 0.5
Perturbation=0.5

printf "Calling calc_sensi.py\n"
python calc_sensi.py $nElements $Perturbation $StartDate $EndDate $JacobianRunsDir $RunName $SensiDir; wait
printf "DONE -- calc_sensi.py\n\n"

#=======================================================================
# Setup GC data directory in workdir
#=======================================================================
GCsourcepth="${JacobianRunsDir}/${RunName}_0000/OutputDir"

printf "Calling setup_GCdatadir.py\n"
python setup_GCdatadir.py $StartDate $EndDate $GCsourcepth $GCDir; wait
printf "DONE -- setup_GCdatadir.py\n\n"

#=======================================================================
# Generate Jacobian matrix files 
#=======================================================================

printf "Calling jacobian.py\n"
python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $FetchTROPOMI; wait
printf " DONE -- jacobian.py\n\n"

#=======================================================================
# Do inversion
#=======================================================================

# Set input values
LonMin={LON_MIN}
LonMax={LON_MAX}
LatMin={LAT_MIN}
LatMax={LAT_MAX}
Res={RES}
posteriorSF="./inversion_result.nc"

printf "Calling invert.py\n"
python invert.py $nElements $JacobianDir $posteriorSF $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $PriorError $ObsError $Gamma $Res; wait
printf "DONE -- invert.py\n\n"

#=======================================================================
# Create gridded posterior scaling factor netcdf file
#=======================================================================
GriddedPosterior="./gridded_posterior.nc"

printf "Calling make_gridded_posterior.py\n"
python make_gridded_posterior.py $posteriorSF $StateVectorFile $GriddedPosterior; wait
printf "DONE -- make_gridded_posterior.py\n\n"

exit 0
