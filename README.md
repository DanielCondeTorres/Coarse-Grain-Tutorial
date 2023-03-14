#  Training school from eCOST: Coarse-Grained Tutorial
This is a short tutorial on how to perform a peptide-membrane interaction simulation using the [Martini](http://www.cgmartini.nl/) force field, together with the [Gromacs](https://www.gromacs.org/) package.

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/199940424-2ad7347e-bbcb-4425-bfc6-bb7a3fca7413.gif" alt="animated"  />
</p>

## Some Theory

### What are computer simulations?

- It is a tool to investigate many particle systems
- It is a quantitave science consist of interplay between experiment and theory
- Useful, because only a few very simplified problems can be solved analytically


### Use computer simulations

- Test models to understand and interpret experimental results
- Test approximate analytical theories
- Explore new regimes hard to reach in experiment
- With a validated model, simulations can be used as 'microscopes' to unveil microscopic mechanisms responsible for certain phenomena

### Types of computer simulations

- [**Monte Carlo simulations:**](https://www.youtube.com/watch?v=WJjDr67frtM) use the knowledge of the hamiltonian and
specific algorithms to appropriately sample the equilibrium microstates
(positions and velocities) of a modeled system. From that sampling we
can obtain equilibrium thermodynamic magnitudes to be compared with
experiment.
- **Molecular Dynamics simulations:** solves the equations of motion of all
the particles in the system to give its temporal evolution given some initial
conditions. **a)** It provides the trajectory (position and velocity as a function of time) of
all the particles in the system. From those, observable quantities can be
computed to compare with experiment. **b)** Can be used to investigate equilibrium (ergodic assumption) and nonequilibrium
phenomena: *The time average of a magnitude equals the average of that magnitude over the [microstates of the same macrostate](https://es.wikipedia.org/wiki/Macroestado_y_microestado)*

### Limitations of classical Molecular Dynamics simulations

- Validity of classical forces and equations if typical distance between particles  is bigger than the De Broglie wavelength
- Validity of intermolecular interactions?
- Limited timescales reachable by current computational capabilities
- Limited size of system by current computational capabilities
- Ergodicity: is the simulated system exploring a relevant fraction of phase space?

### [Periodic boundary conditions](https://www.compchems.com/molecular-dynamics-periodic-boundary-conditions-pbc/) (Born, von Karman (1912), PBC)

To mimic realistic conditions (experiment), in simulations we usually want to describe bulk properties of substances; however, 
due to limited computational power, there is limit in number of particles we can simulate; if walls, a large proportion will be affected by their presence

- No walls, reduces surface effects
- Number density conserved
- Size of simulated box L large compared with potential range
- Suppress long wavelength fluctuations
- Topology of a torus


#### Properties
- Atoms that leave the simulation region at one boundary reenter the region through the opposite face
- Atoms interact with real and image atoms

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212854665-5beaa3aa-bf63-4aa2-9240-f9996f497ffe.gif" alt="animated"  />
</p>

### What is a forcefield?
It is nothing more than a set of parameters and equations that defines the interactions between the particles of our system, based on Newton's second law

## Getting started: Simulation
The first step is to create our simulation box, which will have the following study system:

- Membrane simulating: bacteria (90% POPG and 10% POPE).
- Magainin-2 ETP

The membrane can be created directly with the free software [Charmm-gui](https://www.charmm-gui.org/), all you need to do is register. However, to save time I have moved the work forward and we will be able to download it directly from this repository at the following link:

Magainin-2 will obtain its atomistic representation from the [Protein Data Bank](https://www.rcsb.org/structure/2MAG), if someone is lost, you can get the pdb directly from this repository.
### Load modules that are needed to run and visualize the simulation


In order to use [Gromacs package](https://www.gromacs.org/) 

```
compute -c 16 --mem 32 --gpu
```
```
export OMP_NUM_THREADS=4
```
```
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5
```
and in order to visualize the systems we will use [VMD](https://www.ks.uiuc.edu/Research/vmd/):

```
module load vmd
```
### Atomistic to Coarse-Grain resolution

In this case the representation to be used will be Coarse-Grained, where instead of representing all the atoms in an explicit way we group them in a structure called 'bead' that groups four of these atoms. This way we reduce the number of particles in the system and we can perform computational simulations of larger systems because the computational cost is lower.


<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/213305970-9747cea5-a7ca-485e-a203-1564fb6b0b6c.PNG" />
</p>

To begin with we will use the [Martini] (http://www.cgmartini.nl/) force field, as it allows us to perform this type of simulations, our membrane is already in CG resolution, to transform the peptide (atomistic resolution) to CG we can use the python script provided in Martini's own web page: **martinize.py**


```
python2 martinize.py -f 2mag.pdb -ss HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH -x MAG_CG.pdb -o topol.top -merge all
```
The number of **H** should be the same (or higher) than the BB residues of our ETP.

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
- -nt: allows us to model the peptide in its non-zwiterionic form (in case you want the zwiterionic form, simply do not use the zwiterionic flag)
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


#### IMPORTANT NOTE:

In this tutorial, the membrane has been obtained using the CHARMM-GUI software, which provides us with the .pdb of the membrane, the .itp and a system.top (analogous to topol.top) but with the files we need to define our membrane, the way to proceed would be to include the information of a .top (for example the topol.top) inside the other topology file (system.top), thus grouping all the information in a single file, this is what we will do in the section: **Modifying the topology file.**


 
### Solvating the membrane (not used in this tutorial!!!!)
In case we have the membrane without any type of solvent and we want to add it, it can be used:

```
gmx solvate -cs water.gro -cp caja_grande.pdb -p system.top -scale 4 -o caja_grande_solv.pdb
```

Where the option 'scale 4' allows us not to add the water inside the membrane, but some script can be used to eliminate it.

If we need to create this box to solvate:

```
gmx insert-molecules -ci TFE.pdb  -nmol 1200 -box 5 5 5 -o chx_box.gro
```



<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212669312-70e8a736-a19c-47f0-9426-004ce731e0c1.png" />
</p>

### Insert the molecule in the membrane-water box

In this case we want to insert one Magainin in the simulation box (-nmol 1) in the water solution, so we eliminate waters (W) in order to give space to our peptide

```
gmx insert-molecules -f mem_bact_solvate_clean.pdb -ci MAG_CG.pdb -nmol 1 -replace W -o complete_system.pdb
```



```
grep -wc W complete_system.pdb
```


We add the protein, the correct number of waters and the protein.itp file

you can also use this script provided:


```
python number_of_waters.py -f complete_system.pdb
```

#### Modifying the topology file.

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
*W           27951
*Protein_A   1
 ```             

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212669021-89bb1ddd-68c3-499b-b64e-bf4189a1c377.png" />
</p>

 * means things that we need to modify in system.top

### [Add ions](https://www.researchgate.net/post/Why_do_we_need_to_neutralize_the_charge_of_the_system_by_adding_ions_before_we_proceed_to_energy_minimization_step_in_GROMACS)

 [*Since life does not exist at a net charge, we must add ions to our system.*](https://www.researchgate.net/publication/263952975_Quantifying_Artifacts_in_Ewald_Simulations_of_Inhomogeneous_Systems_with_a_Net_Charge) **See references**

```
gmx grompp -f ions.mdp -c complete_system.pdb -p system.top -o ions.tpr
```


```
gmx genion -s ions.tpr -p system.top -pname NA -nname CL -neutral -o complete_system_ions.gro 
```
We are going to replace the water molecules in order to add the [ions](https://www.researchgate.net/post/Must_system_charge_be_neutral_for_performing_molecular_dynamics_simulation)

Select:
```
Group 4: W
```
**W** means Water

### Minimization

*The purpose of minimization, or relaxation, is to find a local energy minimum of the starting structure so that the molecular dynamics simulation does not immediately “blow up” (i.e. the forces on any one atom are not so large that the atoms move an unreasonable distance in a single timestep). This involves standard minimization algorithms such as steepest descent.* **See references**


```
gmx grompp -f minimization.mdp -c complete_system_ions.gro  -p system.top -o minimize.tpr -maxwarn 1
```



```
gmx_mpi mdrun -v -deffnm minimize 
```

### Make index

We need to create these groups in order to run the .mdp files, because it is going to apply some condictions to the membrane, water, ions and the ETP. When you create the pdbs, you haven't got these groups (you have Cl, Na, POPC ...) but not membrane, ions, etc.

```
gmx make_ndx -f minimize.gro -o index.ndx
```

 We are going to create a Membrane index: POPE and POPG are the lipids:
```
2|3
name 17 Membrane
```
And now, the solvent, therefore: Water and ions

```
4|5
name 18 Water_and_ions
```
In order to leave (finish with the index)

```
q
```
### Equilibration

*Ultimately, we usually seek to run a simulation in a particular thermodynamic ensemble (e.g. the NVE or NVT ensemble) at a particular state point (e.g. target energy, temperature, and pressure) and collect data for analysis which is appropriate for those conditions and not biased depending on our starting conditions/configuration. This means that usually we need to invest simulation time in bringing the system to the appropriate state point as well as relaxing away from any artificially induced metastable starting states. In other words, we are usually interested in sampling the most relevant (or most probable) configurations in the equilibrium ensemble of interest. However, if we start in a less-stable configuration a large part of our equilibration may be the relaxation time (this may be very long for biomolecules or systems at phase equilibrium) necessary to reach the more relevant configuration space.*
**See references**


```
gmx grompp -f equilibration.mdp -c minimize.gro -r minimize.gro -p system.top -o equilibrate.tpr -n index.ndx -maxwarn 1
```

```
gmx_mpi mdrun -v -deffnm equilibrate
```
If it did not load, it is due to a problem because of the iteractive session, so try the following command to run the equilibration.

```
gmx mdrun -v -deffnm minimize -nt 16 -ntmpi 1
```

### Production
*Once equilibration is complete, we may begin collecting data for analysis. Typically this phase is called “production”. The main difference between equilibration and production is simply that in the production simulation, we plan to retain and analyze the collected data. Production must always be preceded by equilibration appropriate for the target production ensemble, and production data should never be collected immediately after a change in conditions (such as rescaling a box size, energy minimizing, or suddenly changing the temperature or pressure) except in very specific applications where this is the goal.* **See references**

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




<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212883913-1ea1135e-fead-4bf0-968a-319cd8d7caff.gif" alt="animated"  />
</p>

#### .mdp files

During the creation of the system, we have used at different points files known as **.mdp files**, these allow us for example which type of algorithms we want to use to find the minimum, during minimisation (such as steepest descent), which type of ensemble production to work on, such as the NVE, NVT or NPT, thus choosing the thermostat, barostat, reference temperatures, the integrator to use to solve the equations of motion, the time step (how often we 'integrate' the equations of motion) or the method to calculate the long range electrostatics such as Ewald, the time step (how often we 'integrate' the equations of motion) or the method to calculate the Long range electrostatics such as Ewald Summation, PME or Reaction-Field, add constraints to the motion of our particles, perform Umbrella-Sampling simulations or how often we check-point our simulation.


### Analysis

Once we go into production, we can study quantitatively how our system behaves, in this section we will present some simple analyses that can be performed directly with GROMACS.

First, we can modify the index and prepare it for the analysis.
For the analysis will be also useful have the protein backbone and a unique water molecule.

```
gmx make_ndx -f minimize.gro -n index.ndx -o index.ndx
```
```
a 6029
a BB
q
```

Where BB is the Backbone of the Protein and 6029 is an unique water molecule, they will be used to calculate the number of contacts and the Radial Distribution Function, respectively.

#### Density profiles


In this first analysis, we will calculate the density of the different components of our system: water and ions (or only water), membrane and protein, we can observe that the value of water is compatible with the theoretical one, at the same time we can observe the density of our lipidic membrane (we can estimate the thickness of the membrane over 4 nm), we also highlight that the ETP is in a stable situation (see the video of the production), resting on the surface of the membrane.

```
gmx density -f trajectory_skip100.xtc -s prod.tpr -n index.ndx -center -ng 3 -o density_skip.xvg
```
```
Select Membrane to center the representation
Select Protein
Select Membrane
Select Water_and_ions
```
```
xmgrace -nxy  density_skip.xvg
```


<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/213301714-215cc30a-0e57-456a-80e3-290585309f93.jpeg" />
</p>


#### Number of contacts

Next, the number of contacts exerted by the beads of the ETP backbone with the membrane has been calculated (using a radius to estimate the contact of 0.6 nm, a typical value for GC simulations) and a total value of approximately 55 contacts with the membrane has been found (Magainin-2 has 23 amino acids = 23 BB beads, and it can be observed that there are beads with more than a single contact). 
 
```
gmx mindist -s prod.tpr -f trajectory_skip100.xtc -on num_contacts_skip.xvg -tu ns -d 0.6 -n index.ndx
```
```
Select BB
Select Membrane
```
```
xmgrace -nxy  num_contacts_skip.xvg
```

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/212896881-5be89031-076a-43d9-a50b-07bf0c62ec41.jpeg" />
</p>

#### Radial distribution function
The following is another typical analysis that students of the Master in Atomistic and Computational Modelling of the UB do when they create their own force field and cry when trying to parallelise it, this is the radial distribution function, which to begin with, depending on the behaviour of the graph we can know if our system is a [solid, liquid or gaseous state](https://commons.wikimedia.org/wiki/File:Simulated_Radial_Distribution_Functions_for_Solid,_Liquid,_and_Gaseous_Argon.svg)

```
gmx rdf -f trajectory_skip100.xtc -s prod.tpr -o RDF_water.xvg -n index.ndx -tu ns
```
```
Select 6029
Select W 
Ctrl+D
```
```
xmgrace -nxy  RDF_water.xvg
```



<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/213301608-1d8b46ba-4e89-4612-afa4-844726e93988.jpeg" />
</p>

The molecules that form this first solvation layer define a space that is inaccessible to other molecules and therefore when we go to greater distances from our central molecule, the number of molecules rapidly decreases and the radial distribution function drops. It then increases until it reaches the second solvation layer of our molecule (and which would form part of the first solvation layer of those molecules). Again, we would have a decrease of the distribution function at larger distances and a new peak. Our reference is located in the central molecule and therefore we can consider it as at rest, but the rest of the molecules are constantly moving. The disordering effect of molecular motion is amplified as we move away from the central molecule (the first layer moves,  the second layer moves relative to the first, and so on), and so for sufficiently large distances, the radial distribution function takes on a uniform value. This same molecular motion explains why the distribution function does not cancel out between the first and second solvation layers. It must be considered that there are configurations in which molecules in the first layer may be moving away from the second layer. That is to say, the probability of finding molecules at that distances is not zero due to the exchange of molecules between the solvation layers.



Finally, I want to point out, that this analysis is also useful to study for example the formation of hydrogen bonds of an amino acid with water molecules of the solvent, as I indicate in the following image, although this analysis would be more useful at an atomistic resolution..., but it serves to clarify the analysis!

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/213033232-cd1af596-c265-43f1-b425-5934b77d24b4.PNG" />
</p>





[Gromacs analysis](https://manual.gromacs.org/documentation/2019/reference-manual/analysis.html) allows you to do a wide variety of analyses, however, we recommend (if you have some notions in Python) the [MDAnalysis](https://www.mdanalysis.org/) module and give free rein to your imagination to study your system.



## References to start with MD simulations

[Braun E, Gilmer J, Mayes HB, Mobley DL, Monroe JI, Prasad S, Zuckerman DM. Best Practices for Foundations in Molecular Simulations [Article v1.0]. Living J Comput Mol Sci. 2019;1(1):5957. doi: 10.33011/livecoms.1.1.5957. Epub 2018 Nov 29. PMID: 31788666; PMCID: PMC6884151](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884151/)

[Gromacs Tutorial](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)

## Importance of colors to make plots for everybody

[ColorBrewer](https://colorbrewer2.org/#type=sequential&scheme=GnBu&n=5)

## Acknowledgments

To **eCOST Action CA21145 EURESTOP European Network for diagnosis and treatment of antibiotic-resistant bacterial infections** for counting on us to give this class.

To **Cesga**, for allowing us the use of their facilities to hold this seminar.

This work has received financial support from the Spanish Agencia Estatal de Investigación (AEI) and the European Regional Development Fund - ERDF (RTI2018-098795-A-I00, PID2019-111327GB-I00, PDC2022-133402-I00 and RYC-2016-20335), by the Xunta de Galicia (ED431B 2022/36, ED431F 2020/05, 06_IN606D_2021_2600276, 02_IN606D_2022_2667887, ED481A 2022/290 and Centro singular de investigación de Galicia accreditation 2019-2022, ED431G 2019/03) and the European Union (ERDF) and by the COST Action CA21145.

I would also like to thank Martín Calvelo, for having taught me with so much patience how to perform my first Molecular Dynamics simulations.

## Team

Daniel Conde-Torres, Alejandro Seco, Fabián Suárez-Lestón, Alfonso Cabezón, Alexandre Blanco, Paula Antelo, Ángel Piñeiro y Rebeca García Fandiño 


## Contact
danielconde.torres@usc.es

alejandro.seco.gonzalez@usc.es

## Social

### LinkedIn

[Daniel Conde-Torres](https://www.linkedin.com/in/daniel-conde-torres-4683b521a)

[Alejandro](https://www.linkedin.com/in/alejandro-seco-gonzalez)

### Instagram

@danicondeofisial

@alexsg14





## Esto parte se implementará en un futuro, pero es para el caso atomístico...

## DSSP

To create the pdbs

```
gmx trjconv -s prod_2.tpr  -f conc.xtc -dt 10 -sep -o traj_.pdb
```

```
pydssp  traj* -o output.result
```
If you have other things than aminoacids... you should go to pdbbio.py:

```
import numpy as np

atomnum = {' N  ':0, ' CA ': 1, ' C  ': 2, ' O  ': 3}

def read_pdbtext_no_checking(pdbstring: str):
    lines = pdbstring.split("\n")
    coords, atoms, resid_old = [], None, None
    for l in lines:
        if l.startswith('ATOM'):
            iatom = atomnum.get(l[12:16], None)
            resid = l[21:26]
            if resid != resid_old:
                if atoms is not None:
                    coords.append(atoms)
                atoms, resid_old = [], resid
            if iatom is not None:
                xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                atoms.append(xyz)
    if atoms is not None:
        coords.append(atoms)
    coords = np.array(coords)
    return coords


def read_pdbtext_with_checking(pdbstring: str):
    lines = pdbstring.split("\n");a=pdbstring.strip()
    coords, atoms, resid_old, check = [], None, None, []
    for l in lines:
        if l[17:20]=='ACE' or l[17:20]=='NH2':#Things that we dont want
            print('Dont use: ',l[17:20])       
        elif l.startswith('ATOM'):
            iatom = atomnum.get(l[12:16], None)
            resid = l[21:26]
            if resid != resid_old:
                if atoms is not None:
                    coords.append(atoms)
                    check.append(atom_check)
                atoms, resid_old, atom_check = [], resid, []
            if iatom is not None:
                xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                atoms.append(xyz)
                atom_check.append(iatom)
    if atoms is not None:
        coords.append(atoms)
        check.append(atom_check)
    coords = np.array(coords)
    # check
    assert len(coords.shape) == 3, "Some required atoms [N,CA,C,O] are missing in the input PDB file"
    check = np.array(check)
    assert np.all(check[:,0]==0), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
    assert np.all(check[:,1]==1), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
    assert np.all(check[:,2]==2), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
    assert np.all(check[:,3]==3), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
    # output
    return coords
```

The script to representate this:

```
import numpy as np
import matplotlib.pyplot as plt
archivo=open('output.result')
file=archivo.readlines()
H_count=[]
E_count=[]
tiempo=[]
for lines in file:
    secondary_structure=lines.split()#leemos las lineas de estructura secundaria
    len_ss=len(secondary_structure[0])#longitud de aminoacidos
    H=0;E=0 #En cada linea establecemos el contador de ss igual a 0
    for elemento in secondary_structure[0]:#contamos la cantidad de elemento por tipo
        if elemento=='H':#alpha helix
            H=H+1     
        elif elemento=='E':#lamina beta
            E=E+1
    E_count.append(E)
    H_count.append(H/len_ss*100)
    

    lista=[0,1,2,3,4,5,6,7,8,9]
    #count the time
    convertir_numero=[]
    for elemento in secondary_structure[1]:
    #    print(elemento)
        try:
            if float(elemento) in lista:
     #           print('Numeros: ',elemento)
                convertir_numero.append(elemento)
            else:
                pass
        except ValueError:
            pass
    tiempo.append(int(''.join(map(str, convertir_numero))))
   # print('wwwwwwwwwwwwwww',int(''.join(map(str, convertir_numero))))

#Ordenamos las listas:
# Creamos una lista de tuplas combinando ambas listas
combinadas = list(zip(tiempo, H_count))

# Ordenamos la lista de tuplas utilizando la lista 1 como clave
combinadas.sort(key=lambda x: x[0])
# Extraemos las listas ordenadas resultantes
tiempo_ordenada, H_count_ordenada = zip(*combinadas)


#vemos relación_frames_tiempo
buscar_relacion_tiempo_frame=open('traj_1.pdb')
#pdb=buscar_relacion_tiempo_frame.readline()
buscar_relacion_tiempo_frame.readline()
pdb=buscar_relacion_tiempo_frame.readline()
factor=float(pdb.split()[5])


'''for line in pdb:
    print(line.split()[3])'''
import matplotlib as mpl
print(tiempo_ordenada[-1])
plt.plot(np.array(tiempo_ordenada)*factor/1000,H_count_ordenada)
plt.ylabel(r'%$\mathbf{(\alpha)-helix}$', fontsize=24)
plt.xlabel(r'$\mathbf{Time (ns)}$', fontsize=24)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.show()
plt.close()
```


