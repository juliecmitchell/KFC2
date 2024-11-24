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
	res_hp=`echo "$LINE"|awk -F"," '{print $1}'`
	pos_per=`echo "$LINE"|awk -F"," '{print $2}'`
	delta_tot=`echo "$LINE"|awk -F"," '{print $3}'`
	core_rim=`echo "$LINE"|awk -F"," '{print $4}'`
	rot5=`echo "$LINE"|awk -F"," '{print $5}'`
	plast4=`echo "$LINE"|awk -F"," '{print $6}'`
	plast5=`echo "$LINE"|awk -F"," '{print $7}'`
	fp10=`echo "$LINE"|sed s/.$// |awk -F"," '{print $8}'`
	echo "-1  1:$res_hp 2:$pos_per 3:$delta_tot 4:$core_rim 5:$rot5 6:$plast4 7:$plast5 8:$fp10" >>$2
done <$1
