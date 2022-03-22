#/bin/bash
#this script was written by Dr. YAN Zeyin at 2020
#contact: zeyin.yan@outlook.com，Southern University of Science and Technology


if [[ $# -eq 1 ]]
then
    logf=$1
    name=${logf%.*}
    steps=`grep -n 'Optimized Parameters' $logf | wc -l`
    #line_nums=`grep -n 'Optimized Parameters' $logf | awk -F ":" '{print $1}'`
    #echo $steps

    line1=1

    num1=`grep -n 'Standard orientation:' $logf | head -n 1 | awk -F ":" '{print $1}'`
    num2=`sed -n $num1",2000p" $logf | grep -n '\---------------------------' | head -n 2 | tail -n 1 | awk -F ":" '{print $1}'`
    num3=`sed -n $num1",2000p" $logf | grep -n '\---------------------------' | head -n 3 | tail -n 1 | awk -F ":" '{print $1}'`
    atoms=$(($num3-$num2-1))
    #echo $atoms
    for (( i=1; i<=$steps; i++ ))
    do
        line2=`grep -n 'Optimized Parameters' $logf | sed -n $i'p' | awk -F ":" '{print $1}'`
        keyline=`sed -n $line1","$line2"p" $logf | grep -n 'Standard orientation:' | tail -n 1 | awk -F ":" '{print $1}'`
        #echo $line1 $line2
        #echo $keyline
        begin=$(($keyline+5))
        end=$(($keyline+$atoms+4))

        #xyz part
        echo $atoms > $name\_$i.xyz
        echo >> $name\_$i.xyz
        sed -n $begin','$end'p' $logf | awk '{print $2,$4,$5,$6}' >> $name\_$i.xyz
    done


else
	echo -e "\033[32mUsage: $0 scan.log \033[0m"
	exit 101
fi

