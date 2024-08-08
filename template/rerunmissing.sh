#/bin/sh -f
# Move the results with missing files to the folder, missing

#option 1: use given list; 
echo usage: ./rerunmissing.sh

#ls -Srlh */*/*/*jump.tif > listjump.txt & #too many files, not working.
# find ./ -type f -name ??_??_?_?_??_??_jump.tif -exec ls -al {} \; | sort -k 5 -n > listjump.txt &
# find ./ -type f -name ??_??_?_?_??_??_jump.tif > listjump.txt &

mkdir missing

for line in `cat listjump.txt`
do
echo $line
listfile=${line/jump.tif/listused.txt}
diri=`echo $(dirname $line)`
nf=`ls $diri |wc -l`  #total number of files should be 8

flag=1
if  [ $nf -ne 8 ] 
then
echo $diri has $nf files, not 8.
flag=0
fi

if [ ! -f $listfile ]
then
    echo "File does not exist in Bash" $listfile
    flag=0
fi

if [ $flag -eq 0 ] 
then
mv $diri  missing/
#exit
fi

#exit
done
