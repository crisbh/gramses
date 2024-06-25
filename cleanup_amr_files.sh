CWD=$(pwd)
TARGETDIR=$1

cd $TARGETDIR

echo "About to clean the following dir:"
pwd

for job in data_job_*;do
    echo "Found $job"
    for dir in $job/output_000*;do 
        echo "Cleaning up files inside $dir" 
        rm -f $dir/amr*;
    done
done

cd $CWD

echo "Done. Back at $CWD and exiting"
