#!/bin/bash

# /!\ FIXME marks places where a custom path, etc. must be specified

#SBATCH --qos=normal
        # Specifies that your job will run in the queue (Quality Of
        # Service) normal, which, if this option is not specified, is
        # the default.
 
#SBATCH --job-name=request-ens
        # Assigns the specified name to the request
 
#SBATCH --output=stdout/request-ens.%j.%N.stdout.log
        # Specifies the name and location of STDOUT where %N is the
        # node where the job runs and %j the job-id. The file will be 
        # written in the workdir directory if it is a relative path.
        # If not given, the default is slurm-%j.out in the workdir.
 
#SBATCH --error=stderr/request-ens.%j.%N.stderr.out
        # Specifies the name and location of STDERR where %N is the
        # node where the job runs and %j the job-id. The file will be 
        # written in the workdir directory if it is a relative path. 
        # If not specified, STDERR will be written to the output file 
        # defined above, or otherwise to slurm-%j.out in the workdir. 
 
#SBATCH --workdir=FIXME
        # Sets the working directory of the batch script before it is
        # executed.
 
#SBATCH --mail-type=END
#SBATCH --mail-user=FIXME
        # Specifies that an email should be sent in case the job fails.
        # Other options include BEGIN, END, REQUEUE and ALL (any state
        # change).

# Take time and date from command line
DATE=$1
TIME=$2
if [ -z "$DATE" ] || [ -z "$TIME" ]
then
	>&2 echo "Must specify arguments YYYY-MM-DD HH."
	exit 1
fi

# Output directory FIXME
TRGD="FIXME/IFS-ENS/ENS-${DATE}T${TIME}Z"
# Output filename prefix
TRGF="${TRGD}/ENS-${DATE}T${TIME}Z-"
# Request file name
RQST="request-ens-${DATE}T${TIME}Z"

mkdir "${TRGD}"

# Wind data on pressure level 400 is not available for older forecasts.
# To avoid the MARS request failing, don't include it in the level list.
# Forecasts earlier than middle of 2008 are missing even more levels
# (I only know that February 2008 has less than September 2008).
if (( "${DATE:0:4}" >= 2010 )); then
	LVLS="10/50/100/200/250/300/400/500/700/850/925/1000"
elif (( "${DATE:0:4}" == 2009 || ( "${DATE:0:4}" == 2008 && ( "${DATE:5:1}" == 1 || "${DATE:6:1}" >= 6 ) ) )); then
	LVLS="10/50/100/200/250/300/500/700/850/925/1000"
else
    LVLS="200/250/300/500/700/850/925/1000"
fi


# 130.128  Temperature
# 131      U-component of wind
# 132      V-component of wind

cat > "${RQST}" <<EOF
retrieve,
    class=od,
    stream=enfo,
    expver=1,
    type=pf,
    grid=1.0/1.0,
    number=1/TO/50,
    step=0/TO/288/BY/6,
    area=90/-180/0/180,
    date=${DATE},
    time=${TIME}:00:00,
    param=130.128,
    levtype=pl,
    levelist=${LVLS},
    target="${TRGF}t.grb"
retrieve,
    param=131,
    levtype=pl,
    levelist=${LVLS},
    target="${TRGF}u.grb"
retrieve,
    param=132,
    levtype=pl,
    levelist=${LVLS},
    target="${TRGF}v.grb"
EOF

exitstatus=0

mars "${RQST}"

if [ $? != 0 ]
then
    echo " The MARS request failed."
    echo
    exitstatus=1
fi

# Convert files to NetCDF and remove original GRIB files
GRIBS="${TRGD}/*.grb"
for GRB in ${GRIBS}
do
	grib_to_netcdf -u time -o "${GRB%.grb}.nc" "${GRB}"  && rm -f "${GRB}"
done
 
exit "$exitstatus"

