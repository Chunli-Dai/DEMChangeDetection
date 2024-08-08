
#wget -i urllist2.txt -o downloadlog

gunzip *.gz

for line in `ls *.tar`; do
echo $line 
tar -xvf $line
ls ${line:0:end-5}*
rm $line
done

