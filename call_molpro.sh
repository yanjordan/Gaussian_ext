#/bin/bash
#this script is to call Molpro WF in DFT function for high level computation in ONIOM
#transport the atomic coordinates to xyz file  called by Molpro input file

if [[ $# == 8 ]]
then
    mpinput=$1
    ncpu=$2
    layer=$3
    input=$4
    output=$5
    msgf=$6
    #cp $input .

    #transfer to Molpro xyz
    gau2mpxyz $input gaump.xyz 

    #lanch molpro computation
    # add the command with a output file
    /share/apps/molpro/molpro_2019_2_linux_x86_64_i8/bin/molpro -n 1 -t $2 --mppx -d$TMPDIR -I$PBS_O_WORKDIR -W$PBS_O_WORKDIR --nobackup --no-xml-output < $1 > molpro.out 

    #get the results from the output file
    energy=`grep "EMBED/" molpro.out | awk -F "=" '{print $2}'`
    printf "%20.12e%20.12e%20.12e%20.12e\n" $energy 0.0 0.0 0.0 > $output
    #numatom=9
    line=`grep -n "PROJ EMBED" molpro.out | tail -n 1 | awk -F ":" '{print $1}'`
    line1=`grep -n "Nuclear force" molpro.out  | tail -n 1 | awk -F ":" '{print $1}'`
    begin=$(($line+4))
    end=$(($line1-2))
    #begin=$(($line+4))
    #end=$(($line+3+$numatom))
    #sed -n $begin','$end'p' molpro.out | awk '{printf "%20.12e%20.12e%20.12e" $2 $3 $4}' >> gauout.dat
    grad=`sed -n $begin','$end'p' molpro.out | awk '{print $2, $3, $4}'`
    printf "%20.12e%20.12e%20.12e\n" $grad >> $output

    sed -i "s/e/D/g" $output

else
    echo -e 'call_molpro.sh input.com ncpus'
    echo -e 'More infor see in Gaussian xx website for more information.'
fi
