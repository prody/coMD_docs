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
  
  #Copy the conf file to these new folder
  cp -avf tmd/4AKE.conf tmd$i/.
  cp -avf tmd/4AKE.job tmd$i/.

  cp -avf tmdr/1AKE.conf tmdr$i/.
  cp -avf tmdr/1AKE.job tmdr$i/.	

  #Generate the current minimization folder
  mkdir tmdm$i
  mkdir tmdmr$i

  #Copy the conf file into these new folder
  cp -avf tmdm/4AKEm.conf tmdm$i/.
  cp -avf tmdm/4AKEm.job tmdm$i/.

  #Copy the conf file into these new folder
  cp -avf tmdmr/1AKEm.conf tmdmr$i/.
  cp -avf tmdmr/1AKEm.job tmdmr$i/.

  echo $i > deger.dat

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

 
  #MATLAB generated structure is copied to common2
  cp -avf main2/final_structure.dcd common2/.
  cp -avf main2r/final_structurer.dcd common2r/.

  #Generate the target file for TMD simualations
  vmd -dispdev text -e  tmd_target_structure.pgn
  vmd -dispdev text -e  tmd_target_structurer.pgn

  wait

  #Copy this structure to the current folder
  cp -avf adjust.pdb tmd$i/.
  cp -avf adjustr.pdb tmdr$i/.

  #Run TMD



cd tmd$i

JOBID=`qsub 4AKE.job | cut -f 3 -d' '`

# http://tldp.org/LDP/abs/html/

echo "The job id is $JOBID"
cd ..

cd tmdr$i

JOBID2=`qsub 1AKE.job | cut -f 3 -d' '`

# http://tldp.org/LDP/abs/html/

echo "The second job id is $JOBID2"
cd ..


while qstat -j $JOBID > /dev/null 2>&1
do
  sleep 60
  echo "Job still running"
done
echo "Job complete"

 
while qstat -j $JOBID2 > /dev/null 2>&1
do
  sleep 60
  echo "Second job still running"
done
echo "Second job complete"

 

  cp -avf tmd$i/res.coor tmd/.
  cp -avf tmd$i/res.vel tmd/.
  cp -avf tmd$i/res.xsc tmd/.

  cp -avf tmdr$i/res.coor tmdr/.
  cp -avf tmdr$i/res.vel tmdr/.
  cp -avf tmdr$i/res.xsc tmdr/.

cd tmdm$i

JOBID=`qsub 4AKEm.job | cut -f 3 -d' '`

# http://tldp.org/LDP/abs/html/

echo "The job id is $JOBID"

cd ..

cd tmdmr$i

JOBID2=`qsub 1AKEm.job | cut -f 3 -d' '`

# http://tldp.org/LDP/abs/html/

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
  echo "Second job still running"
done
echo "Second job complete"


  #Copy the final structure to the common folder
  cp -avf tmdm$i/res.coor common2/.
  cp -avf tmdm$i/res.vel common2/.
  cp -avf tmdm$i/res.xsc common2/.


  cp -avf tmdmr$i/res.coor common2r/.
  cp -avf tmdmr$i/res.vel common2r/.
  cp -avf tmdmr$i/res.xsc common2r/.


#Check convergence

vmd -dispdev text -e  RMSD_checker.pgn

wait


#Stop coMD simulation of convergence is achieved

if (( `cat status.dat` ))
then
  exit
fi
  
done




#sleep 5

killall vmd_LINUXAMD64

echo "Script Bitti"
