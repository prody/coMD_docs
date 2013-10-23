#!/bin/bash

clear
echo " coMD Started"


#Use conformers in 500i and 500f as END structures A and B

 cp -avf 500i/res.coor common2/.
 cp -avf 500i/res.vel common2/.
 cp -avf 500i/res.xsc common2/.


 cp -avf 500f/res.coor common2r/.
 cp -avf 500f/res.vel common2r/.
 cp -avf 500f/res.xsc common2r/.

wait

#Compile Matlab

cd main2
  mcc -m main2
wait
cd ..

cd main2r
  mcc -m main2r
wait
cd ..


# Cycles of coMD 

for i in {1..100}


do

  #Perform alignments with targets
  vmd -dispdev text -e  align_start_with_final_structure.pgn
  vmd -dispdev text -e  align_start_with_final_structurer.pgn
  wait

  #Generate pdb files with the CA coordinates of the targets only
  vmd -dispdev text -e  CA_final_structure.pgn
  vmd -dispdev text -e  CA_final_structurer.pgn
  wait

  echo "Running for i: $i"
  #Generate the current folder TMD folder
  mkdir tmd$i
  mkdir tmdr$i

  cp -avf target.pdb tmd$i/.
  cp -avf targetr.pdb tmdr$i/
  
	
  #Generate the current minimization folder
  mkdir tmdm$i
  mkdir tmdmr$i

  

  cd main2

  JOBID=`qsub main2.job | cut -f 3 -d' '`

  # http://tldp.org/LDP/abs/html/

  echo "The job id is $JOBID"
  cd ..

  cd main2r
 
  JOBID2=`qsub main2r.job | cut -f 3 -d' '`

  echo "The second job id is $JOBID2"
  cd ..

  while qstat -j $JOBID > /dev/null 2>&1
  do
  sleep 60
  echo "Job Still running"
  done
  echo "Job complete"

  while qstat -j $JOBID2 > /dev/null 2>&1
  do
  sleep 60
  echo "Second job Still running"
  done
  echo "Second job complete"



  # Copy all files required for TMD simulations
  cp -avf main2/cg.dcd tmd$i/.
  cp -avf main2/final_structure.dcd tmd$i/.
  cp -avf main2/oran.dat tmd$i/.
  cp -avf main2/main2.out tmd$i/.

  # Copy all files required for TMD simulations
  cp -avf main2r/cgr.dcd tmdr$i/.
  cp -avf main2r/final_structurer.dcd tmdr$i/.
  cp -avf main2r/oranr.dat tmdr$i/.
  cp -avf main2r/main2r.out tmdr$i/.

 
  
#Stop coMD simulation of convergence is achieved

if (( `cat status.dat` ))
then
	if (( `cat statusr.dat` ))
	then	
  		exit
	fi
fi
  
done




#sleep 5

killall vmd_LINUXAMD64

echo "Script Bitti"
