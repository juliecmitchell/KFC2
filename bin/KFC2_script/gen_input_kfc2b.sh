#!/bin/sh
# used to generate train data for SVM and for k-con
usage()
{
	echo "usage: $0 input.txt output"
}
if [ $# -ne 2 ]
then
	usage
	exit 1
fi
while read LINE
do
	res_size=`echo "$LINE"|awk -F"," '{print $1}'`
	pos_per=`echo "$LINE"|awk -F"," '{print $2}'`
	ratio5=`echo "$LINE"|awk -F"," '{print $3}'`
	rot4=`echo "$LINE"|awk -F"," '{print $4}'`
	rot5=`echo "$LINE"|awk -F"," '{print $5}'`
	hp5=`echo "$LINE"|awk -F"," '{print $6}'`
	fp9=`echo "$LINE"|sed s/.$// |awk -F"," '{print $7}'`
	echo "-1  1:$res_size 2:$pos_per 3:$ratio5 4:$rot4 5:$rot5 6:$hp5 7:$fp9" >>$2
done <$1
