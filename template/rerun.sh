#/bin/sh -f

#grep 'DUE TO TIME LIMIT' *err  > df; mv df log_rerun.txt

#Automate the following steps
# grep -B 1 60778779 outrun1
# sbatch qsubj188.pbs

input="log_rerun.txt"
file1="outrun1"
ns=1 #start number, default 1

#for line in `cat log_rerun.txt` #read space by space
#do
#done

# read line by line
i=0
while IFS= read -r line
do
let i=i+1

#echo $line
jobid=${line:0:8}
str=`grep -B 1 $jobid outrun1 |grep joblist`
str1=`echo ${str:end-15:15}`
jobno=`echo $str1 |awk -Fjoblist '{print $2}'`
jobname=qsubj$jobno.pbs
echo echo $i jobid jobname $jobid $jobname
if [ $i -ge $ns ] 
then
  echo hi $i
  sbatch $jobname
fi

done < "$input"

