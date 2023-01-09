# Coarse-Grain-Tutorial
This is a short tutorial on how to perform a peptide-membrane interaction simulation using the Martini force field, together with the Gromacs package.

<p align="center">
  <img src="https://user-images.githubusercontent.com/117435891/199940424-2ad7347e-bbcb-4425-bfc6-bb7a3fca7413.gif" alt="animated"  />
</p>

# Getting started
The first step is to create our simulation box, which will have the following study system:

- Membrane simulating that of a healthy mammal: (100% POPC).
- Magainin-2 ETP

The membrane can be created directly with the free software [Charmm-gui](https://www.charmm-gui.org/), all you need to do is register. However, to save time I have moved the work forward and we will be able to download it directly from this repository at the following link:

Magainin-2 will obtain its atomistic representation from the [Protein Data Bank](https://www.rcsb.org/structure/2MAG), if someone is lost, you can get the pdb directly from this repository.
# Load modules that are needed to run and visualize the simulation
In order to use Gromacs package 
'''
module load gcc/system openmpi/4.0.5 gromacs/2021-PLUMED-2.7.1 mdanalysis/1.1.1
'''
and in order to visualize the systems:
'''
module load vmd
'''
# Atomistic to Coarse-Grain resolution
