#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 4000
#SBATCH -t 0-03:00
#SBATCH -o logfiles/inversion_script_%j.out
#SBATCH -e logfiles/inversion_script_%j.err

#=======================================================================
# Configuration
#=======================================================================
NCLUST={CLUSTERS}
STARTDAY={START}
ENDDAY={END}
JACRUNSDIR={RUNDIRS}
SENSIDIR="./inversion/Sensi"
GCDATADIR="./inversion/data_GC"
JACOBIANDIR="./inversion/data_converted"
CLUSTERSPTH="{MYPATH}/input_data_permian/Clusters_permian_kmeans.nc"

# This doesn't matter until we want to run posterior simulations after the inversion
POSTRUNDIR="${JACRUNSDIR}/posterior_run"

# Load the Python environment
# Comment this out for now and instruct users to activate their own Python envs
# module load Anaconda3/5.0.1-fasrc01
# source activate tropomi_oilgas_jupyter_3.6

# Make sure specified paths exist
if [[ ! -d ${JACRUNSDIR} ]]; then
    echo "\n${JACRUNSDIR} does not exist. Please fix JACRUNDIR."
    exit 1
fi
if [[ ! -f ${JCLUSTERSPTH} ]]; then
    echo "\n${CLUSTERSPTH} does not exist. Please fix CLUSTERSPTH.."
    exit 1
fi

#=======================================================================
# Postprocess the SpeciesConc and LevelEdgeDiags files from GEOS-Chem
#=======================================================================

# Only matters for Kalman filter multi-week inversions
firstsimswitch="True"

# Irrelevant unless running posterior simulation after inversion
postdir="${POSTRUNDIR}/{RUNNAME}_Posterior_0000"

python postproc_diags.py $JACRUNSDIR $postdir $STARTDAY $firstsimswitch; wait
echo "Postprocessed SpeciesConc and LevelEdgeDiags files, FSS=$firstsimswitch"

#=======================================================================
# Calculate GEOS-Chem sensitivities and save to Sensi directory
#=======================================================================

python calc_sensi.py $STARTDAY $ENDDAY $JACRUNSDIR $SENSIDIR; wait
echo "Saved GC sensitivity files to $SENSIDIR"

#=======================================================================
# Setup GC data directory in workdir
#=======================================================================
GCsourcepth="${JACRUNSDIR}/{RUNNAME}_0000/OutputDir"

python setup_GCdatadir.py $STARTDAY $ENDDAY $GCsourcepth $GCDATADIR; wait
echo "Setup GEOS-Chem data directory for hourly data files"

#=======================================================================
# Generate Jacobian matrix files 
#=======================================================================
jacstoredir="${JACOBIANDIR}"
mkdir $jacstoredir

python step1_jacobian.py $STARTDAY $ENDDAY $jacstoredir; wait
echo "Generated Jacobian matrix files in $jacstoredir"

#=======================================================================
# Do inversion
#=======================================================================
posteriorSF_outputpath="./inversion_result.nc"

python step2_invert.py $posteriorSF_outputpath $jacstoredir; wait
echo "Did inversion and saved results to $posteriorSF_outputpath"

#=======================================================================
# Create gridded posterior scaling factor netcdf file
#=======================================================================
gridded_posterior_path="./gridded_posterior.nc"

python make_gridded_posterior.py $posteriorSF_outputpath $CLUSTERSPTH $gridded_posterior_path; wait
echo "Mapped posterior scaling factors to domain grid"

exit 0
