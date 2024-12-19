from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.constraints import FixAtoms
from ase.constraints import FixedLine
from ase.md import MDLogger
from ase.io import read, write
from ase.io.extxyz import read_extxyz
import numpy as np
import time
import os
import ssl
from mace.calculators.mace import MACECalculator


############################
# UPDATE OUTPUTS
############################
# Define the path to the file
file_path = 'RESTART_COUNT'

def read_restart_count():                                                                                                                                                                                                             
    # Open the file to read the last entry
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Extract the last entry if the file is not empty
        last_entry = float(lines[-1].strip()) if lines else 0
    print(f"RESTART_COUNT initial value: {last_entry}")
    return last_entry

# Execute the function to update the RESTART_COUNT
init_restart = read_restart_count()


############################
# SET SIMULATION TIMES
############################
equilibration_1_steps = 45000
equilibration_2_steps = 45000
production_steps = 400000
timestep = 0.5
production_track = (production_steps * timestep)/1000

############################
# EQUILIBRATION 1
############################
# Initial configuration
init_conf = read('AAw12.2-17-1.0-water-material-1h3o.pdb', '0')
init_conf.set_calculator(MACECalculator('/mnt/lustre/e1000/home/ec214/ec214/xr223_cirrus/models/c-h2o-def/c-h2o-def_swa.model', device='cuda', default_dtype='float32'))

# Freeze all C  positions
fixed_atoms = FixAtoms(mask=[atom.symbol == 'C' for atom in init_conf])
init_conf.set_constraint(fixed_atoms)

# Equilibration at 300K
MaxwellBoltzmannDistribution(init_conf, temperature_K=300)
dyn = Langevin(init_conf, 0.5 * units.fs, temperature_K=300, friction=0.0025 / units.fs)
def write_frame():
    dyn.atoms.write('equi-1-traj.xyz', append=True)
dyn.attach(write_frame, interval=200)
dyn.attach(MDLogger(dyn, init_conf, 'equi-1-md.log', header=True, stress=False, peratom=False, mode="a"), interval=200)
dyn.run(equilibration_1_steps)

write('end.pdb', init_conf)

############################
# EQUILIBRATION 2
############################
# Initial configuration
init_conf =  read('equi-1-traj.xyz', '-1')
init_conf.set_calculator(MACECalculator('/mnt/lustre/e1000/home/ec214/ec214/xr223_cirrus/models/c-h2o-def/c-h2o-def_swa.model', device='cuda', default_dtype='float32'))

velocities = init_conf.get_velocities()
init_conf.set_velocities(velocities)

# Freeze one C atom L1, freeze xy one C atom L2
fixed_bl_xyz = FixAtoms([atom.index for atom in init_conf if ((atom.symbol == 'C') and (atom.index == 0))])
fixed_ul_xy = FixedLine([atom.index for atom in init_conf if ((atom.symbol == 'C') and (atom.index == 4))], direction=[0, 0, 1])
init_conf.set_constraint([fixed_bl_xyz, fixed_ul_xy])

# Equilibration at 300K
#MaxwellBoltzmannDistribution(init_conf, temperature_K=300)
dyn = Langevin(init_conf, 0.5 * units.fs, temperature_K=300, friction=0.0025 / units.fs)
def write_frame():
    dyn.atoms.write('equi-2-traj.xyz', append=True)
dyn.attach(write_frame, interval=200)
dyn.attach(MDLogger(dyn, init_conf, 'equi-2-md.log', header=True, stress=False, peratom=False, mode="a"), interval=200)
dyn.run(equilibration_2_steps)

write('end-2.pdb', init_conf)

############################
# PRODUCTION
############################
# Initial configuration
init_conf =  read('equi-2-traj.xyz', '-1')
init_conf.set_calculator(MACECalculator('/mnt/lustre/e1000/home/ec214/ec214/xr223_cirrus/models/c-h2o-def/c-h2o-def_swa.model', device='cuda', default_dtype='float32'))

velocities = init_conf.get_velocities() 
init_conf.set_velocities(velocities)

# Freeze one C atom L1, freeze xy one C atom L2
fixed_bl_xyz = FixAtoms([atom.index for atom in init_conf if ((atom.symbol == 'C') and (atom.index == 0))])
fixed_ul_xy = FixedLine([atom.index for atom in init_conf if ((atom.symbol == 'C') and (atom.index == 4))], direction=[0, 0, 1])
init_conf.set_constraint([fixed_bl_xyz, fixed_ul_xy])

# Production run at 300K
#MaxwellBoltzmannDistribution(init_conf, temperature_K=300)
dyn = Langevin(init_conf, 0.5 * units.fs, temperature_K=300, friction=0.0025 / units.fs)
def write_frame():
    dyn.atoms.write(f'prod-traj-{init_restart}.xyz', append=True)
dyn.attach(write_frame, interval=20)
dyn.attach(MDLogger(dyn, init_conf, f'prod-md-{init_restart}.log', header=True, stress=False, peratom=False, mode="a"), interval=20)
dyn.run(production_steps)

write(f'end-prod-{init_restart}.pdb', init_conf)



############################
# UPDATE OUTPUTS
############################
# Define the path to the file
file_path = 'RESTART_COUNT'

def update_restart_count(production_track):
    try:
        # Open the file to read the last entry
        with open(file_path, 'r') as file:
            lines = file.readlines()
            # Extract the last entry if the file is not empty
            last_entry = float(lines[-1].strip()) if lines else 0
    except FileNotFoundError:
        # If the file doesn't exist, start with 0
        last_entry = 0

    # Add 500 to the last entry
    new_entry = last_entry + production_track

    # Append the new entry to the file
    with open(file_path, 'a') as file:
        file.write(f"{new_entry}\n")

    print(f"Updated RESTART_COUNT with new value: {new_entry}")

# Execute the function to update the RESTART_COUNT
update_restart_count(production_track)

hint_file='SWITCH_RES'
with open(hint_file, 'a') as file:
        file.write(f"File printed to let HPC know we should restart") 
