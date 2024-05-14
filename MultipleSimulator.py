import os
import numpy as np
import pandas as pd
import time
import pathlib

HOME = str(pathlib.Path.home())

# The directory of gromacs
gmx22 = f"{HOME}/Gromacs2022.3/Installation/bin/gmx"

# Directing to the Working and Saving Directory(put your .top and .gro and .cpt files here)
print(os.getcwd())
print("Before running the simulations, you should do the below instructions:\n\n"
      "\t1) Copy the .gro, .top, and .cpt file in the Test2 folder or the WORKING directory.(especially files from the NPT process)\n"
      "\t2) Delete or move everything from the WORKING and SAVING directories.\n"
      "\t(Test2 folder is for your input and output, and Test3 foleder is for your text terminal saving)\n\n")
dir1 = input("Where is your WORKING and INPUT and OUTPUT directory?(ex: ~/Documents/Test2 == 'def')\n->")
if dir1 == "def":
    dir1 = f"{HOME}/Documents/Test2"
dir2 = input("Where is your SAVING directory?(ex: ~/Documents/Test3 == 'def')\n->")
if dir2 == "def":
    dir2 = f"{HOME}/Documents/Test3"


# The function to append and writing for Saving directory
def Write(name, content, aow):
    """

    :param name: name of the file with extention(text.txt)
    :param content: the content that wrote in the file
    :param aow: you can choose a or w to append or create a file
    :return: saving a content to a 'name' file
    """
    f = open(name, aow)
    f.write(content)
    f.close()


# Making the Different Folders
name1 = input("What is your PDB protein name?\t\t\t->")           # The Protein name for naming our folders
name2 = input("What time is it?(ex: 00.12.30)\t\t\t->")           # The Date for naming our folders
name3 = input("What is the Beginning time in nanoseconds?\t->")     # For determining that if we have in-continuous simulation or the beginning simulation.
                                                                        # And for naming our folders.
name4 = input("What is the ending time in nanoseconds?\t\t->")      # For naming our folders.
name5 = input("Which number is the file?\t\t\t->") 		        # The number in front them.


# The Accelerations numbers you want
acc_list = []
acc = input("Which accelerations do you want?\n(type 'all' if you want the following Accelerations:\n"
            "['0.00', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '0.30']\n"
            "But if you want some specific Accelerations value to be in simulation,\n"
            "type 'num' for add number by number and press Enter and Enter the numbers individually.\n"
            "But if you want to add a whole numbers as list, type 'list' and write numbers and seperate them with space. If you want to end it, then you can press Enter without any number to type to end it.)\n->")
if acc == "all":
    acc_list = ['0.00', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '0.30']
elif acc == "num":
    while acc != "":
        acc_list.append(acc)
        if acc_list[-1] == "num":
            acc_list.pop(-1)
        acc = input("Enter your number and press Enter,\n our you don't have any more number,\n"
                    "press just Enter.\n->")
    print(f"Accelerations: {acc_list}")
elif acc == "list":
    acc = input("Enter your list by typing the number and SPACE it to differentiate the numbers.\nlike this:\n0.00 1.00 1.44 5.11\n")
    index_list = [-1]
    for i in range(len(acc)):
        if acc[i] == " ":
            acc_list.append(acc[index_list[-1]+1 : i])
            index_list.append(i)
    acc_list.append(acc[index_list[-1] + 1:])
    print(f"Accelerations: {acc_list}")

# First time snap
tFirst = time.asctime()

# Making the essential root for all the simulations and folders
os.chdir(dir1)
os.system(f"mkdir {name1}_{name3}to{name4}ns_{name5}")
# Making the folders after this root
dir3 = f"{dir1}/{name1}_{name3}to{name4}ns_{name5}"
os.chdir(dir3)
print(os.getcwd())
for i in range(len(acc_list)):
    os.system(f'mkdir -p "({acc_list[i]}){name2}_{name1}"')


# Making the .mdp files in each and copy the .top and .cpt and .gro file in each.
for i in range(len(acc_list)):
    os.chdir(f'{dir3}/({acc_list[i]}){name2}_{name1}')

    # Making the .mdp file
    print(f"md.acc_{acc_list[i]}.mdp  ----> The mdp files\n\n")        # The Helping Part
    time.sleep(1)
    f = open(f"md.acc_{acc_list[i]}.mdp", "w")
    f.write(f"title                   = {name1} NPT equilibration \n"
            "; Run parameters\n"
            "integrator              = md        ; leap-frog integrator\n"
            f"nsteps                  = {int((float(name4)-float(name3))*500000)}   ; 2 * {int((float(name4)-float(name3))*500000)} = {int((float(name4)-float(name3))*1000)} ps ({int((float(name4)-float(name3)))} ns)\n"
            "dt                      = 0.002     ; 2 fs\n"
            "; Output control\n"
            "nstxout                 = 0         ; suppress bulky .trr file by specifying \n"
            "nstvout                 = 0         ; 0 for output frequency of nstxout,\n"
            "nstfout                 = 0         ; nstvout, and nstfout\n"
            "nstenergy               = 5000      ; save energies every 10.0 ps\n"
            "nstlog                  = 5000      ; update log file every 10.0 ps\n"
            "nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps\n"
            "compressed-x-grps       = System    ; save the whole system\n"
            "; Bond parameters\n"
            "continuation            = yes       ; Restarting after NPT \n"
            "constraint_algorithm    = lincs     ; holonomic constraints \n"
            "constraints             = h-bonds   ; bonds involving H are constrained\n"
            "lincs_iter              = 1         ; accuracy of LINCS\n"
            "lincs_order             = 4         ; also related to accuracy\n"
            "; Neighborsearching\n"
            "cutoff-scheme           = Verlet    ; Buffered neighbor searching\n"
            "ns_type                 = grid      ; search neighboring grid cells\n"
            "nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme\n"
            "rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)\n"
            "rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)\n"
            "; Electrostatics\n"
            "coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics\n"
            "pme_order               = 4         ; cubic interpolation\n"
            "fourierspacing          = 0.16      ; grid spacing for FFT\n"
            "; Temperature coupling is on\n"
            "tcoupl                  = V-rescale             ; modified Berendsen thermostat\n"
            "tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate\n"
            "tau_t                   = 0.1     0.1           ; time constant, in ps\n"
            "ref_t                   = 300     300           ; reference temperature, one for each group, in K\n"
            "; Pressure coupling is on\n"
            "pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT\n"
            "pcoupltype              = isotropic             ; uniform scaling of box vectors\n"
            "tau_p                   = 2.0                   ; time constant, in ps\n"
            "ref_p                   = 1.0                   ; reference pressure, in bar\n"
            "compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1\n"
            "; Periodic boundary conditions\n"
            "pbc                     = xyz       ; 3-D PBC\n"
            "; Dispersion correction\n"
            "DispCorr                = EnerPres  ; account for cut-off vdW scheme\n"
            "; Velocity generation\n"
            "gen_vel                 = no        ; Velocity generation is off \n"
            "; Non-equilibrium MD stuff\n"
            "acc-grps                 = SOL\n"
            f"accelerate               = {float(acc_list[i])} 0.0 0.0\n"
            "cos-acceleration         = 0\n")
    f.close()

    # Coping the .gro, .top and .cpt files here
    print(f"{dir3}/({acc_list[i]}){name2}_{name1}  ---> The WORKING directory\n\n")                      # The Helping Part
    time.sleep(1)
    os.system(f'cp {dir1}/*.gro "{dir3}/({acc_list[i]}){name2}_{name1}"')
    os.system(f'cp {dir1}/*.top "{dir3}/({acc_list[i]}){name2}_{name1}"')
    os.system(f'cp {dir1}/*.cpt "{dir3}/({acc_list[i]}){name2}_{name1}"')


# The Simulations and the noting(Saving) process
print(f"{dir2}/{name2}_{name1}_{name5}in   ----> The Main folder name\n\n")      # The Helping Part
time.sleep(1)
Write(f"{dir2}/{name2}_{name1}_{name5}in", "", "w")

for i in range(len(acc_list)):
    print(f"{dir2}/{name2}_{name1}_{name5}in   ----> The Main folder name\n\n")  # The Helping Part
    time.sleep(1)
    Write(f"{dir2}/{name2}_{name1}_{name5}in",
          "-------------------------------------------------------------------------------------------------------------------\n"
          f"1) md.acc_{acc_list[i]}.mdp:\n"
          "\n"
          f"title                   = {name1} NPT equilibration \n"
          "; Run parameters\n"
          "integrator              = md        ; leap-frog integrator\n"
            f"nsteps                  = {int((float(name4)-float(name3))*500000)}   ; 2 * {int((float(name4)-float(name3))*500000)} = {int((float(name4)-float(name3))*1000)} ps ({int((float(name4)-float(name3)))} ns)\n"
          "dt                      = 0.002     ; 2 fs\n"
          "; Output control\n"
          "nstxout                 = 0         ; suppress bulky .trr file by specifying \n"
          "nstvout                 = 0         ; 0 for output frequency of nstxout,\n"
          "nstfout                 = 0         ; nstvout, and nstfout\n"
          "nstenergy               = 5000      ; save energies every 10.0 ps\n"
          "nstlog                  = 5000      ; update log file every 10.0 ps\n"
          "nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps\n"
          "compressed-x-grps       = System    ; save the whole system\n"
          "; Bond parameters\n"
          "continuation            = yes       ; Restarting after NPT \n"
          "constraint_algorithm    = lincs     ; holonomic constraints \n"
          "constraints             = h-bonds   ; bonds involving H are constrained\n"
          "lincs_iter              = 1         ; accuracy of LINCS\n"
          "lincs_order             = 4         ; also related to accuracy\n"
          "; Neighborsearching\n"
          "cutoff-scheme           = Verlet    ; Buffered neighbor searching\n"
          "ns_type                 = grid      ; search neighboring grid cells\n"
          "nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme\n"
          "rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)\n"
          "rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)\n"
          "; Electrostatics\n"
          "coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics\n"
          "pme_order               = 4         ; cubic interpolation\n"
          "fourierspacing          = 0.16      ; grid spacing for FFT\n"
          "; Temperature coupling is on\n"
          "tcoupl                  = V-rescale             ; modified Berendsen thermostat\n"
          "tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate\n"
          "tau_t                   = 0.1     0.1           ; time constant, in ps\n"
          "ref_t                   = 300     300           ; reference temperature, one for each group, in K\n"
          "; Pressure coupling is on\n"
          "pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT\n"
          "pcoupltype              = isotropic             ; uniform scaling of box vectors\n"
          "tau_p                   = 2.0                   ; time constant, in ps\n"
          "ref_p                   = 1.0                   ; reference pressure, in bar\n"
          "compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1\n"
          "; Periodic boundary conditions\n"
          "pbc                     = xyz       ; 3-D PBC\n"
          "; Dispersion correction\n"
          "DispCorr                = EnerPres  ; account for cut-off vdW scheme\n"
          "; Velocity generation\n"
          "gen_vel                 = no        ; Velocity generation is off \n"
          "; Non-equilibrium MD stuff\n"
          "acc-grps                 = SOL\n"
          f"accelerate               = {float(acc_list[i])} 0.0 0.0\n"
          "cos-acceleration         = 0\n"
          "\n"
          "2) The Grompp command: \n"
          f"$ gmx22 grompp -f md.acc_{acc_list[i]}.mdp -c npt.gro -t npt.cpt -p topol.top -o md.acc_{acc_list[i]}.tpr\n\n"
          
          "-In the Terminal:\n", "a")

    print(f"{dir3}/({acc_list[i]}){name2}_{name1}   ----> The each Folder in the Main Folder\n\n")  # The Helping Part
    time.sleep(1)
    os.chdir(f"{dir3}/({acc_list[i]}){name2}_{name1}")
    print(f"{gmx22} grompp -f md.acc_{acc_list[i]}.mdp -c *.gro -t *.cpt -p topol.top -o md.acc_{acc_list[i]}.tpr"
              f" 2>&1 | tee -a {dir2}/{name2}_{name1}_{name5}in   ----> The grompp command\n\n")  # The Helping Part
    time.sleep(1)
    os.system(f"{gmx22} grompp -f md.acc_{acc_list[i]}.mdp -c *.gro -t *.cpt -p topol.top -o md.acc_{acc_list[i]}.tpr"
              f" 2>&1 | tee -a {dir2}/{name2}_{name1}_{name5}in")
    print(f"{dir2}/{name2}_{name1}_{name5}in  -------> The SAVING directory\n\n")   # The Helping Part
    time.sleep(1)
    Write(f"{dir2}/{name2}_{name1}_{name5}in",
          "-The files that are made:\n"
          f"md.acc_{acc_list[i]}.tpr\n\n"
          "3) The Mdrun command:\n"
          f"$ gmx22 mdrun -v -deffnm md.acc_{acc_list[i]}\n"
          "-The Terminal:\n\n", "a")
    print(f"{gmx22} mdrun -v -deffnm md.acc_{acc_list[i]} 2>&1 | tee -a {dir2}/{name2}_{name1}_{name5}in  ---> the mdrun command\n\n")    # The Helping Part
    time.sleep(1)
    os.system(f"{gmx22} mdrun -v -deffnm md.acc_{acc_list[i]} 2>&1 | tee -a {dir2}/{name2}_{name1}_{name5}in")
    Write(f"{dir2}/{name2}_{name1}_{name5}in",
          f"The files that made:\n"
          f"md.acc_{acc_list[i]}.cpt\n"
          f"md.acc_{acc_list[i]}.edr\n"
          f"md.acc_{acc_list[i]}.gro\n"
          f"md.acc_{acc_list[i]}.log\n"
          f"md.acc_{acc_list[i]}.xtc\n\n\n", "a")

#########################################################
# The SAVING file correction!!!!

os.chdir(dir2)
print(os.getcwd())

text = open(f"{name2}_{name1}_{name5}in", "r")  # Reading the input SAVING file
texts = text.readlines()
text.close()

print("acc>>>>>", acc_list)
ans1 = len(acc_list)    # number of simulation for number of iteration for this instruction.
print(ans1)

textout = open(f"{name2}_{name1}_{name5}", "w+")  # The new output SAVING file

text4find = []
for i in range(len(texts)):   # Making a SAVING file for reading the lines from 0 to 11 word
    try:
        text4find.append(texts[i][:11])
    except IndexError:
        text4find.append(texts[i])

fin = -1
lisbeg = []
lisfin = []
for j in range(ans1):                          # The 'step 1000' and 'step 499900' lines and send them to a list for finding indeces
    print(int(float(name3)))
    print(fin + 1)
    beg = text4find.index(f'step {(int(float(name3))) + 200}, w', fin+1)
    print(beg)
    lisbeg.append(beg)
    print(int(float(name4)))
    print(fin + 1)
    fin = text4find.index(f'step {(int(float(name4)) * 500000) - 100}', fin+1)
    print(fin)
    lisfin.append(fin)

lisbeg.append(len(texts)-1)
print(lisbeg)
print(lisfin)

i = -1
while i < len(lisfin):                  # The Correction part
    print(i, len(lisfin))
    if i == -1:
        for line in range(0, lisbeg[0]):
            textout.write(texts[line])
        i += 1
    else:

        for line in range(lisfin[i], lisbeg[i+1]):
            textout.write(texts[line])
        i += 1


textout.close()
# Final time snap
tFinal = time.asctime()

print("------------------THE TIME------------------")
print(f"Beginning time:\t\t {tFirst}\n"
      f"Final time:\t\t {tFinal}")

