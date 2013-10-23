.. CoMD documentation master file, created by
   sphinx-quickstart on Wed Oct  9 15:02:08 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


#################################
|CM| Tutorial
#################################

=================================
What you need to start
=================================

**Having NAMD, VMD and MATLAB installed**

* NAMD needs to be installed on your cluster/supercomputer to execute the Targeted Molecular Dynamics (TMD) and Energy Minimization Simulations.
* MATLAB needs to be installed on your cluster/supercomputer  to execute the MC/Metropolis Scheme.
* VMD needs to be installed on your cluster/supercomputer  for the communication between MATLAB and NAMD.

Please see the schematic figure (`Figure 2`_ in the `Paper`_)  below for a visual description where each program is required.

.. _Figure 2: http://www.cell.com/biophysj/image/S0006-3495%2813%2900934-X?imageId=gr2&imageType=large
.. _Paper: http://www.cell.com/biophysj/fulltext/S0006-3495%2813%2900934-X

**Having Equilibrated Conformers of your End States**

You took your end states, structure A and B, from the Protein Data Bank. What you need to do now is to equilibrate them in Molecular Dynamics (MD), preferably in the presence of explicit solvent (and also membrane if a membrane protein is simulated).

To be consistent with the TMD+MD scheme of the |CM| cycles, I prefer to perform a 500 step minimization prior to generating the first set of ANM modes.

.. |CM| replace:: \ *co*\MD


.. figure:: Pictures/Webpage_protocol.jpg
   :scale: 25 %
   
   Figure 1: Schematic Description of |CM|
   	
   	The programs that are used to excecute each part are shown in the figure


======================================
Download the Script
======================================

Download the :download:`initializing code <./Initialize/Initialize.tar.gz>` and the main :download:`|CM| code <./MainCode/coMD.tar.gz>`.

Take a quick look into these folders. 

What subfolders do you have? What type of files are in those folders? Obtain an overview of the files.


======================================
 Prepare your files
======================================

As mentioned before you need to equilibrate your states **A** and **B**. Once you have done that minimize them for 500 steps and put the restart files into the folders **500i** and **500f**.

======================================
 What is the MC/Metropolis Parameter **a**, and how is it adjusted
======================================

Please open the file main2/main2.m inside the |CM| folder. 

In order to attain a specific MC/Metropolis acceptance ratio the parameter **a** inside the MC/Metropolis condition (see below) has to be set::

	exp(-(En-Ep)*a)>rand()

This is performed by first initialzing **a**.  

**a** is adjusted inside the code as follows::

	if prevper(2)>0.95
    		a=prevper(5)*1.5;
	elseif prevper(2)<0.85
    		a=prevper(5)/1.5;
	end
   

Here 0.85 and 0.95 refer to the lower and upper limits of a %90 MC/Metropolis Acceptance ratio simulations.

It the acceptance ratio becomes larger than %95 or lower than %85 then **a** is readjusted.

======================================
 Initialize the MC/Metropolis Parameters
======================================


To initialize  **a** we perform |CM| simulation without the *TMD+MD* part.

The script used to initialize the parameters for AK can be found :download:`here  <./Initialize/Initialize.tar.gz>`.

Go into this folder and simply submit the shell script main.sh by typing the following into the terminal::

	nohup sh ./main.sh > main.out &


Once the initialization simulations has converged take the **oran.dat** and **oranr.dat** files inside the folders main2 and main2r, respectively, and put them into your main  |CM| simulation folder (under main2 and main2r).

======================================
 Take some time to check what pgn files do
======================================
 
There are four types of pgn files that coMD uses::

	align_start_with_final_structure.pgn

	CA_final_structure.pgn

	tmd_target_structure.pgn

	RMSDchecker.pgn


*align_start_with_final_structure.pgn* aligns the current conformer with its target.

*CA_final_structure.pgn* generates the target pdb's for MC/Metropolis Simulations with CA coordinates only.

*tmd_target_structure.pgn* generates the target pdb's for Targeted Molecular Dynamics Simulations. 

PLEASE NOTE THAT YOU NEED PDB FILES (CA.pdb and CAr.pdb) THAT CONTAINS AS MANY CA ATOMS AS YOUR MC/METROPOLIS ALGORITHM USES.

Above the pgn files for the A-->B direction were discussed. There are the same files with "r" endings for the B-->A direction.


======================================
 Number of residues
======================================
Please note the number of CA atoms has to be the same for states A and B.

If they are not, then you need to change the pgn script (See below) so that the selections have the same number of CA atoms.

This is really much easier than it sounds. You simply define instead of::

	 [atomselect  "name CA"]

the following::

	 [atomselect "resid xxx to yyy and name CA"]



======================================
 You are all set to run the |CM| simulation
======================================

Go into the |CM| folder and simply submit the shell script main.sh by typing the follwing into the terminal::

	nohup sh ./main.sh > main.out &


======================================
 pgn files to generate transition trajectories from intermediates
======================================

:download:`transformation.pgn  <./Investigate/transformation.pgn>` and :download:`transformationr.pgn  <./Investigate/transformationr.pgn>` can be used to generate dcd files 
for the A-->B and B-->A trajectories.

Please note that intermediate conformers in those dcd files are ordered such that the first frame of the B-->A trajectory
follows the last frame of the A-->B trajectory.

In addition you have to adjust these pgn files with respect to the total number of |CM| cycles. It is described inside the pgn files.

.. toctree::
   :maxdepth: 2

   credits



