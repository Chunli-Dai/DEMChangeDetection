#/bin/sh -f
# 
# Modified from /home/dai.56/chunli/scripts/run_overlap_strip.sh 
# Modified based on run_change2.sh (November 2020)
# Changes: run 25 subtiles in one job to reduce the total number of jobs. 

inputtype=3; 
inputtype=1; 
inputtype=5;
# 1 %based on input xid etc.
# 2 %find the block based on input coordinates;
# 3 %Use a rectangle box around the input coordinates;
# 5 %use a list of tile names

if [ $inputtype -eq 2 -o $inputtype -eq 3  ] ; then
#loneq=;lateq=;workdir=site1;
file='aoi.txt'
nline=`wc -l < $file`
echo Total number of sites: $nline.
#5 11 12 16 17 18 22 24 25
for (( i=11; i<=$nline; i++ )) #1:80 $nline
#for i in 1 22 25
do 
latloneq=`sed -n ''$i'p' aoi.txt |awk  -F' ' '{ print $1, $2 }'`

workdir=site$i
if [ ! -d $workdir ] 
then
    mkdir $workdir
fi
cd $workdir
ln -fs ../mat0.mat .
ln -fs ../*.m .
ln -fs ../job*pbs .

echo $inputtype > input.txt
echo $latloneq >> input.txt
jobid=`qsub -N job1 jobpar.pbs`

#wait until this job is done to run the second job
#while false
while true
do
sleep 5s #wait 5 seconds
out=`qstat $jobid`
status=`echo $out|awk '{print $19}'`
#echo $status 
if [[ "$status" == "C" ]]
then
   break
fi
done
echo $jobid $status "finished; start to work on next site."

#pid="$!"
#wait $pid
#tail --pid=$pid -f /dev/null
cd ../
#exit
done #i

#else #inputtype=1 every 2km size grid.
elif [ $inputtype -eq 1 ] ; then
#40_18_2_2_2m_v3.0 40_18_2_1_2m_v3.0
#40_18_1_2_2m_v3.0 40_18_1_1_2m_v3.0 39_18_1_2_2m_v3.0
#40_17_2_2_2m_v3.0 40_17_2_1_2m_v3.0 39_17_2_2_2m_v3.0
#40_17_1_2_2m_v3.0 40_17_1_1_2m_v3.0 39_17_1_2_2m_v3.0

#for xid=1:22;
#for yid=37:74; #Alaska
# arcticdem_07_canada_ellesmere: yid 27:35; xid:26:39

for (( xid=26; xid<=39; xid++ )) #1:80; 33:34
do 
for (( yid=27; yid<=35; yid++ )) #1:80 32:33
do 
for (( xids=1; xids<=2; xids++ )) #1:2
do 
for (( yids=1; yids<=2; yids++ ))  #1:2
do 

xidc=`printf "%2.2d" "$xid"`
yidc=`printf "%2.2d" "$yid"`
xidsc=`printf "%1.1d" "$xids"`
yidsc=`printf "%1.1d" "$yids"`

workdirtile1="$yidc"_"$xidc"_"$xidsc"_"$yidsc"
echo workdirtile1 $workdirtile1
if [ ! -d $workdirtile1 ] 
then
    mkdir $workdirtile1
fi

cd $workdirtile1
pwd
coderundir=`pwd`; # Code with updated parameters for a specific job.
pwdsv=`pwd`;

jobidg=(); # reset $jobid matrix
for (( xidss=1; xidss<=25; xidss++ )) #1:25
do 

xidssc="`printf "%2.2d" "$xidss"`"
str1="$yidc"_"$xidc"_"$xidsc"_"$yidsc"_"$xidssc" #group name # e.g. 41_16_1_1_11
# run the following 25 subtiles sequentially in one job.
mkdir $str1
cd $str1
pwd
ln -fs ../../job_group.pbs .

#tmpstr=`qsub -F "$str1" ./job_group.pbs`
#tmpstr=`sbatch --export=str1="$str1" ./job_group.slurm `
#jobidg[$xidss-1]=`echo $tmpstr|awk '{print $4}'`
#jobidg[$xidss-1]=`qsub -F "$str1" ./job_group.pbs `
#example: qsub -v str1='41_16_1_1_11' ./job_group.pbs
jobidg[$xidss-1]=`qsub ./job_group.pbs `
#exit

cd ../

done #xidss
cd ../
#exit #32_33_2_2

echo list all jobs ${jobidg[@]} #list all jobs

#jobid=`qsub -N job1 jobpar.pbs`
jobid=${jobidg[24]}
echo The last job ID: $jobid.

#wait until the last job is done to run the second job
#while false
while true
do
sleep 5s #wait 5 seconds
out=`qstat $jobid`
#out=`squeue --job $jobid`
#echo out $out
status=`echo $out|awk '{print $19}'`
#status=`echo $out|awk '{print $13}'`
#echo $status
if [ "$status" == "C" -o -z "$out" ]
then
   echo job $jobid finished.
   break
fi
done #while true
echo $jobid $status "finished; start to work on next tile."

done #yids
done #xids
done #yid
done  #xid

elif [ $inputtype -eq 4 ] ; then
#;%use boundary
# see /home/dai.56/chunliwork/ice/greenland/run_change2.sh
echo $inputtype


elif [ $inputtype -eq 5 ] ; then
#use a list of tile names
#tilelist=('39_17_2_2' '40_17_1_2' '40_17_1_1' '39_17_1_2');
tilelist=('39_17_1_2');

# get tilelist from a filelist
# ls -d ??_??_?_? > tilelist
i=0
for line in `cat tilelist`
do
i=$i+1;
tilelist[$i-1]=$line;
done

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

shrundir=`pwd`; #e.g. /u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/
rm joblist[0-9]*
count=0
countj=1;

for (( i=1; i<=$nlist; i++ ))
do
tilename=${tilelist[$i-1]} #39_17_2_2 "$yidc"_"$xidc"_"$xidsc"_"$yidsc"
xidc=${tilename:3:2}
yidc=${tilename:0:2}
#xidsc=${tilename:6:1}
#yidsc=${tilename:8:1}

for (( xids=1; xids<=2; xids++ )) #1:2
do 
for (( yids=1; yids<=2; yids++ ))  #1:2
do 

xidsc=`printf "%1.1d" "$xids"`
yidsc=`printf "%1.1d" "$yids"`

workdirtile1="$yidc"_"$xidc"_"$xidsc"_"$yidsc"
echo workdirtile1 $workdirtile1
if [ ! -d $workdirtile1 ] 
then
    mkdir $workdirtile1
fi

# 39_17_2_2
cd $workdirtile1
pwd
#pwdsv=`pwd`;
dir_tile50km=`pwd`;

jobidg=(); # reset $jobid matrix
for (( xidss=1; xidss<=25; xidss++ )) #1:25
do 

xidssc="`printf "%2.2d" "$xidss"`"
str1="$yidc"_"$xidc"_"$xidsc"_"$yidsc"_"$xidssc" #group name # e.g. 41_16_1_1_11
# run the following 25 subtiles sequentially in one job.
unlink $str1
mkdir $str1
# 41_16_1_1_11
cd $str1
pwd
ln -fs $shrundir/job_group.pbs .

#jobidg[$xidss-1]=`qsub ./job_group.pbs `
#exit #hi

#bundle every 40 jobs.
let "count+=1"
echo $count
if [[ $count -lt 40 ]] ; then
ofile=$shrundir/joblist$countj # /u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/joblist1
echo $ofile
strtmp=`pwd`
echo "$strtmp"'/ ' >> "$ofile"

else #==40 #when we collect 40 jobs, we submit them
ofile=$shrundir/joblist$countj
strtmp=`pwd`
echo "$strtmp"'/ ' >> $ofile

#submit job
pbsfile=$shrundir/qsubj$countj.pbs #/u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/qsubj1.pbs
cp /home/chunlidai/blue/apps/landslide/template/parallel.sbatch $pbsfile
newtext=$ofile
oldtext="/home/chunli/chunliwork/work/landslide/testparallel/joblist"
echo $oldtext $newtext  
#use # as separator in sed.
sed -i 's#'$oldtext'#'$newtext'#g' $pbsfile
lastsubtdir=`pwd` #last subtile directory of a job list.
cd $shrundir  #so that output .err .out would be in the main work directory.
sbatch $pbsfile 
cd $lastsubtdir #compatible with old logic.
#reset
let "countj+=1"
count=0

#exit  #hi

fi #40


#wait if the current number of jobs in queue is >=450
while true
do
sleep 1s #wait 5 seconds
njobs=`squeue -u chunlidai | wc -l`
#maxjobs=450;
maxjobs=14; #16; # max: 9*40 for 5 years; 18 *40 until may 2023
let "njobs-=1"
if [[ $njobs -lt $maxjobs ]]
then
   echo The number of current jobs is $njobs ', < '$maxjobs . 
   break
else
   echo The number of current jobs is $njobs ', wait to submit new jobs until it is <' $maxjobs .
   sleep 10m #wait 1 minute
fi
done
# wait

cd $dir_tile50km/ #../  #41_16_1_1/

done #xidss
cd $shrundir/ #../  #$shrundir

echo list all jobs ${jobidg[@]} #list all jobs

done #yids
done #xids

done # tile list for (( i=1; i<=$nlist; i++ ))

fi

#submit the last job
pbsfile=$shrundir/qsubj$countj.pbs #/u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/qsubj1.pbs
cp /home/chunlidai/blue/apps/landslide/template/parallel.sbatch $pbsfile
newtext=$ofile
oldtext="/home/chunli/chunliwork/work/landslide/testparallel/joblist"
echo $oldtext $newtext  
#use # as separator in sed.
sed -i 's#'$oldtext'#'$newtext'#g' $pbsfile
lastsubtdir=`pwd` #last subtile directory of a job list.
cd $shrundir  #so that output .err .out would be in the main work directory.
sbatch $pbsfile
cd $lastsubtdir #compatible with old logic.

