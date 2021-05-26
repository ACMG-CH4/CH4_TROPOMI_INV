#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 4000
#SBATCH -t 0-03:00
#SBATCH -o run_inversion_%j.out
#SBATCH -e run_inversion_%j.err

#=======================================================================
# Configuration
#=======================================================================
NCLUST={CLUSTERS}
STARTDAY={START}
ENDDAY={END}
MYPATH={MY_PATH}
RUNNAME={RUN_NAME}
SPINUPDIR="${MYPATH}/${RUNNAME}/spinup_run"
JACRUNSDIR="${MYPATH}/${RUNNAME}/jacobian_runs"
POSTRUNDIR="${MYPATH}/${RUNNAME}/posterior_run"
CLUSTERSPTH="${MYPATH}/input_data_permian/Clusters_permian_kmeans.nc"
SENSIDIR="./Sensi"
GCDATADIR="./data_GC"
JACOBIANDIR="./data_converted"

# Only matters for Kalman filter multi-week inversions
firstsimswitch=true

# Load the Python environment
# Comment this out for now and instruct users to activate their own Python envs
# module load Anaconda3/5.0.1-fasrc01
# source activate tropomi_oilgas_jupyter_3.6

# Make sure specified paths exist
if [[ ! -d ${JACRUNSDIR} ]]; then
    echo "${JACRUNSDIR} does not exist. Please fix JACRUNDIR."
    exit 1
fi
if [[ ! -f ${CLUSTERSPTH} ]]; then
    echo "${CLUSTERSPTH} does not exist. Please fix CLUSTERSPTH."
    exit 1
fi

#=======================================================================
# Postprocess the SpeciesConc and LevelEdgeDiags files from GEOS-Chem
#=======================================================================

echo "Calling postproc_diags.py, FSS=$firstsimswitch"
if "$firstsimswitch"; then
    if [[ ! -d ${SPINUPDIR} ]]; then
	echo "${SPINUPDIR} does not exist. Please fix SPINUPDIR or set firstsimswitch to False."
	exit 1
    fi
    PREVDIR=$SPINUPDIR
else
    PREVDIR=%$POSTRUNDIR
    if [[ ! -d ${POSTRUNDIR} ]]; then
	echo "${POSTRUNDIR} does not exist. Please fix POSTRUNDIR."
	exit 1
    fi
fi
echo "  - Hour 0 for ${STARTDAY} will be obtained from ${PREVDIR}"
    
python postproc_diags.py $RUNNAME $JACRUNSDIR $PREVDIR $STARTDAY; wait
echo "DONE -- postproc_diags.py"
echo ""

#=======================================================================
# Calculate GEOS-Chem sensitivities and save to Sensi directory
#=======================================================================

# Temporary input definition
PERTURBATION=0.5

echo "Calling calc_sensi.py"
python calc_sensi.py $NCLUST $PERTURBATION $STARTDAY $ENDDAY $JACRUNSDIR $RUNNAME $SENSIDIR; wait
echo "DONE -- calc_sensi.py"
echo ""

#=======================================================================
# Setup GC data directory in workdir
#=======================================================================
GCsourcepth="${JACRUNSDIR}/${RUNNAME}_0000/OutputDir"

echo "Calling setup_GCdatadir.py"
python setup_GCdatadir.py $STARTDAY $ENDDAY $GCsourcepth $GCDATADIR; wait
echo "DONE -- setup_GCdatadir.py"
echo ""

#=======================================================================
# Generate Jacobian matrix files 
#=======================================================================

echo "Calling jacobian.py"
python jacobian.py $STARTDAY $ENDDAY; wait
echo " DONE -- jacobian.py"
echo ""

#=======================================================================
# Do inversion
#=======================================================================

# New inputs (for Permian)
LON_MIN=-111
LON_MAX=-95
LAT_MIN=25
LAT_MAX=39
PRIOR_ERR=0.5
OBS_ERR=15
GAMMA=0.25

echo "Calling invert.py"
python invert.py $NCLUST $JACOBIANDIR $LON_MIN $LON_MAX $LAT_MIN $LAT_MAX $PRIOR_ERR $OBS_ERR $GAMMA; wait
echo "DONE -- invert.py"
echo ""

#=======================================================================
# Create gridded posterior scaling factor netcdf file
#=======================================================================
posteriorSF="./inversion_result.nc"
gridded_posterior="./gridded_posterior.nc"

echo "Calling make_gridded_posterior.py"
python make_gridded_posterior.py $posteriorSF $CLUSTERSPTH $gridded_posterior; wait
echo "DONE -- make_gridded_posterior.py"
echo ""

exit 0
