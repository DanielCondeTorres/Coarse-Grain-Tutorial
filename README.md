


# Coarse-Grain-Tutorial
This is a short tutorial on how to perform a peptide-membrane interaction simulation using the Martini force field, together with the Gromacs package.

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/199940424-2ad7347e-bbcb-4425-bfc6-bb7a3fca7413.gif" alt="animated"  />
</p>


# Getting started
The first step is to create our simulation box, which will have the following study system:

- Membrane simulating: bacteria (90% POPG and 10% POPE).
- Magainin-2 ETP

The membrane can be created directly with the free software [Charmm-gui](https://www.charmm-gui.org/), all you need to do is register. However, to save time I have moved the work forward and we will be able to download it directly from this repository at the following link:

Magainin-2 will obtain its atomistic representation from the [Protein Data Bank](https://www.rcsb.org/structure/2MAG), if someone is lost, you can get the pdb directly from this repository.
# Load modules that are needed to run and visualize the simulation
In order to use [Gromacs package](https://www.gromacs.org/) 

```
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5
```
and in order to visualize the systems we will use [VMD](https://www.ks.uiuc.edu/Research/vmd/):

```
module load vmd
```
# Atomistic to Coarse-Grain resolution

In this case the representation to be used will be Coarse-Grained, where instead of representing all the atoms in an explicit way we group them in a structure called 'bead' that groups four of these atoms. This way we reduce the number of particles in the system and we can perform computational simulations of larger systems because the computational cost is lower.

To begin with we will use the [Martini] (http://www.cgmartini.nl/force) field, as it allows us to perform this type of simulations, our membrane is already in CG resolution, to transform the peptide (atomistic resolution) to CG we can use the python script provided in Martini's own web page: **martinize.py**


```
python2 martinize.py -f 2mag.pdb -ss HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH -x MAG_CG.pdb -o topol.top -merge all
```
- -ss: indicates the secondary structure that we want, in this casa α-helix
1.  "F": "Collagenous Fiber",                                                                 
2.  "E": "Extended structure (beta sheet)",                                                    
3.  "H": "Helix structure",                                                                    
4.  "1": "Helix start (H-bond donor)",                                                         
5.  "2": "Helix end (H-bond acceptor)",                                                        
6.  "3": "Ambivalent helix type (short helices)",                                              
7.  "T": "Turn",                                                                               
8.  "S": "Bend",                                                                               
9.  "C": "Coil", 
- -x: is the CG structure
- -merge: Merge chains
- -topol.top: defines our system, the files that we are going to use and also the molecules of our system:
 ```
#include "martini.itp"
#include "Protein_A.itp"
[ system ]
; name
Martini system from 2mag(3).pdb
[ molecules ]
; name        number
Protein_A        1
 ```
 1. "martini.itp": parameters and equations that describes our forcefield
 2. "Protein_A.itp": angles, distances, dihedrals... parameters that define our particle 
 
# Solvating the membrane (not used in this tutorial)
In case we have the membrane without any type of solvent and we want to add it, it can be used:

```
gmx solvate -cs water.gro -cp caja_grande.pdb -p system.top -scale 4 -o caja_grande_solv.pdb
```

Where the option 'scale 4' allows us not to add the water inside the membrane, but some script can be used to eliminate it.





<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212669312-70e8a736-a19c-47f0-9426-004ce731e0c1.png" />
</p>

# Insert the molecule in the membrane-water box

In this case we want to insert one Magainin in the simulation box (-nmol 1) in the water solution, so we eliminate waters (W) in order to give space to our peptide

```
gmx insert-molecules -f mem_bact_solvate_clean.pdb -ci MAG_CG.pdb -nmol 1 -replace W -o complete_system.pdb
```

```
python number_of_waters.py -f complete_system.pdb
```
We add the protein, the correc number of waters and the protein itp file





```
#include "martini_v2.2.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_ions.itp"
*#include "Protein_A.itp"

[ system ]
; name
Martini system in water

[ molecules ]
; name        number
POPE 25
POPG 225
POPE 25
POPG 225
*W           27966
*Protein_A   1
 ```             

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212669021-89bb1ddd-68c3-499b-b64e-bf4189a1c377.png" />
</p>

# Add ions

 *Since life does not exist at a net charge, we must add ions to our system.* **See references**

```
gmx grompp -f ions.mdp -c complete_system.pdb -p system.top -o ions.tpr
```


```
gmx genion -s ions.tpr -p system.top -pname NA -nname CL -neutral -o complete_system_ions.gro 
```
We are going to replace the water molecules in order to add the ions

Select:
```
Group 4: W
```
**W** meand Water

# Minimization

*The purpose of minimization, or relaxation, is to find a local energy minimum of the starting structure so that the molecular dynamics simulation does not immediately “blow up” (i.e. the forces on any one atom are not so large that the atoms move an unreasonable distance in a single timestep). This involves standard minimization algorithms such as steepest descent.* **See references**


```
gmx grompp -f minimization.mdp -c complete_system_ions.gro  -p system.top -o minimize.tpr -maxwarn 1
```



```
gmx mdrun -v -deffnm minimize
```

# Make index

We need to create these groups in order to run the .mdp files, because it is going to apply some condictions to the membrane, water, ions and the ETP. When you create the pdbs, you haven't got these groups (you have Cl, Na, POPC ...) but not membrane, ions, etc.

```
gmx make_ndx -f minimize.gro -o index.ndx
```

POPE and POPG

```
2|3
name 17 Membrane
4|5
name 18 Water_and_ions
```

# Equilibration

*Ultimately, we usually seek to run a simulation in a particular thermodynamic ensemble (e.g. the NVE or NVT ensemble) at a particular state point (e.g. target energy, temperature, and pressure) and collect data for analysis which is appropriate for those conditions and not biased depending on our starting conditions/configuration. This means that usually we need to invest simulation time in bringing the system to the appropriate state point as well as relaxing away from any artificially induced metastable starting states. In other words, we are usually interested in sampling the most relevant (or most probable) configurations in the equilibrium ensemble of interest. However, if we start in a less-stable configuration a large part of our equilibration may be the relaxation time (this may be very long for biomolecules or systems at phase equilibrium) necessary to reach the more relevant configuration space.*
**See references**


```
gmx grompp -f equilibration.mdp -c minimize.gro -r minimize.gro -p system.top -o equilibrate.tpr -n index -maxwarn 1
```

```
gmx mdrun -v -deffnm equilibrate
```

# Production

```
gmx grompp -f run.mdp -c equilibrate.gro -p system.top -o prod.tpr -n index -maxwarn 1
```

```
#!/bin/bash
#SBATCH -t  00:30:00 # execution time. Ex: 1/2 hour
#SBATCH --mem-per-cpu=1GB
#SBATCH -n 1 -c 1# number of tasks, number of cores
#SBATCH --ntasks-per-node=1
module load gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
srun gmx_mpi mdrun -pin on -cpi -noappend -s prod.tpr                                                                                                     
```


```
sbatch prod.sh
```


# References to start with MD simulations

[Braun E, Gilmer J, Mayes HB, Mobley DL, Monroe JI, Prasad S, Zuckerman DM. Best Practices for Foundations in Molecular Simulations [Article v1.0]. Living J Comput Mol Sci. 2019;1(1):5957. doi: 10.33011/livecoms.1.1.5957. Epub 2018 Nov 29. PMID: 31788666; PMCID: PMC6884151](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884151/)

[Gromacs Tutorial](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)

# Importance of colors to make plots for everybody

[ColorBrewer](https://colorbrewer2.org/#type=sequential&scheme=GnBu&n=5)
