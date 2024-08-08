#/bin/sh -f

#option 1: use given list; 
echo usage: ./removesml.sh

#ls -Srlh */*/*/*jump.tif > listjump.txt & #too many files, not working.
# find ./ -type f -name ??_??_?_?_??_??_jump.tif -exec ls -al {} \; | sort -k 5 -n > listjump.txt &

for line in `cat smljump.txt`
#for line in `cat jump_bad_region18.txt`
do

echo $line
rm ${line:0:end-8}*

done
