print("Requirments:\n"
      " Installing following packages in the non-conda environment:\n"
      "     1) pandas\t\t$ pip install pandas\n"
      "     2) xlsxwriter\t$ pip install xlsxwriter\n"
      "     3) openpyxl\t$ pip install openpyxl\n"
      "     4) pymol\t\t$ conda install -c conda-forge -c schrodinger pymol-bundle\n"
      "     5) mdtraj\t\t$ conda install -c conda-forge mdtraj\n"
      "     6) contact_map\t$ conda install -c conda-forge contact_map\n"
      "                 \t or\n"
      "                 \t$ pip install contact_map\n"
      "     7) matplotlib\t$ pip install matplotlib\n")

import os
import pandas as pd
import matplotlib.pyplot as plt
import xlsxwriter as xw
import time
import pymol as pm
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
from contact_map import OverrideTopologyContactDifference as otcd

# Constants:
nameOfUser = ""
gmxVersion = "2022.3"
gmxExe = f"/home/{nameOfUser}/Gromacs{gmxVersion}/Installation/bin/gmx"

# Functions:
def copyFiles(filenameWithDir, finalDir):
    os.system(f"cp {filenameWithDir} {finalDir}")
def moveFiles(filenameWithDir, finalDir):
    os.system(f"mv {filenameWithDir} {finalDir}")
def listOfFilesAndFolders():
    listOfName = os.popen("ls").readlines()
    listOfOutputs = []
    for name in listOfName:
        listOfOutputs.append(name[:-1])
    return listOfOutputs
def ansi(color, text, background=0, style=49):
    """

    :param color:   31 = red
                    32 = green
                    33 = yellow
                    34 = cobalt blue
                    35 = pink
                    36 = cyan
                    37 = white
                    39 = gray
                    94 = blue navy
    :param text: Enter your text.
    :param background: ...
    :param style:   0 = normal
                    1 = a little brighter
                    2 = a lighter
                    3 = with underline
                    4 = a little lighter
                    5 = the background is the color
    :return:
    """
    return f"\33[{background};{style};{color}m{text}\33[0m"
def The_Processor(n, i):
    """

    :param n: The number of processes.
    :param i: The process that we are in right now.
    :return: The process bar.
    """
    import time

    if i / n == 0:
        print("The Process:")
        print(f"|{ansi(34, '                    ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 0 < (i / n) * 100 <= 5:
        print("The Process:")
        print(f"|{ansi(34, '█                   ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 5 < (i / n) * 100 <= 10:
        print("The Process:")
        print(f"|{ansi(34, '██                  ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 10 < (i / n) * 100 <= 15:
        print("The Process:")
        print(f"|{ansi(34, '███                 ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 15 < (i / n) * 100 <= 20:
        print("The Process:")
        print(f"|{ansi(34, '████                ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 20 < (i / n) * 100 <= 25:
        print("The Process:")
        print(f"|{ansi(34, '█████               ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 25 < (i / n) * 100 <= 30:
        print("The Process:")
        print(f"|{ansi(34, '██████              ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 30 < (i / n) * 100 <= 35:
        print("The Process:")
        print(f"|{ansi(34, '███████             ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 35 < (i / n) * 100 <= 40:
        print("The Process:")
        print(f"|{ansi(34, '████████            ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 40 < (i / n) * 100 <= 45:
        print("The Process:")
        print(f"|{ansi(34, '█████████           ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 45 < (i / n) * 100 <= 50:
        print("The Process:")
        print(f"|{ansi(34, '██████████          ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 50 < (i / n) * 100 <= 55:
        print("The Process:")
        print(f"|{ansi(34, '███████████         ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 55 < (i / n) * 100 <= 60:
        print("The Process:")
        print(f"|{ansi(34, '████████████        ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 60 < (i / n) * 100 <= 65:
        print("The Process:")
        print(f"|{ansi(34, '█████████████       ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 65 < (i / n) * 100 <= 70:
        print("The Process:")
        print(f"|{ansi(34, '██████████████      ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 70 < (i / n) * 100 <= 75:
        print("The Process:")
        print(f"|{ansi(34, '███████████████     ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 75 < (i / n) * 100 <= 80:
        print("The Process:")
        print(f"|{ansi(34, '████████████████    ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 80 < (i / n) * 100 <= 85:
        print("The Process:")
        print(f"|{ansi(34, '█████████████████   ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 85 < (i / n) * 100 <= 90:
        print("The Process:")
        print(f"|{ansi(34, '██████████████████  ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    elif 90 < (i / n) * 100 <= 95:
        print("The Process:")
        print(f"|{ansi(34, '███████████████████ ')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")
    else:
        print("The Process:")
        print(f"|{ansi(34, '████████████████████')}| {ansi(32, f'{int((i / n) * 100)}%')} to complete")

    time.sleep(0.6)


print(f"***********************************\n"
      f"* {ansi(95, 'Welcome to my ANALYSIS program.')} *\n"
      f"***********************************\n")

print(f"\n {ansi(31, 'Please read the instruction before choosing one Analysis method:', style = 7)}\n"
      f"* Each Protein have a folder (the '{ansi(94, 'Protein Folder')}' that contain folders that repeat series of same Simulation (the '{ansi(94, 'Repetition Folder')}')\n"
      f"    Each 'repetition folder' contain a series of folders that each contain a simulation of interested acceleration ('{ansi(94, 'acceleration folder')}')\n"
      f"    Each 'Acceleration Folder' contain .xtc and .tpr files and following data files in it.\n\n"
      f"- Make three folders in 'Documents' path with names: 'Test1', 'Test2', 'Test3'\n"
      f"- Clean all folders and files from Test1, Test2, Test3.\n"
      f"- Make a folder named '{ansi(37, 'proteins', style=3)}' in Test2 folder.\n"
      f"- Put all the 'Repetition Folder's from all 'Protein Folder's in the 'proteins' folder.\n"
      f"- Run the Program.\n"
      f"- You have three address you need to input them:\n"
      f"        .Permenant Output ADDRESS(Test1 or final_Saving_path\n"
      f"        .Temporary Output ADDRESS(Test3 or Saving_path\n"
      f"        .Input files ADDRESSS(Test2 or source_path\n")

alys_name = "something"
while alys_name != "":
    alys_name = input(f"What {ansi(1, 'analysis')} do you want?\n(The {ansi(5, 'options', style = 3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style = 7)} to exit the program) \n")
    print(alys_name)

    # -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
    if alys_name == "dssp":  # The DSSP analysis
        ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'READ ME', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        print(f"\n                          {ansi(32, 'DSSP Introduction', style=4)}\n"
              f" The {ansi(31, '310 helix')}, α {ansi(31, 'helix')} and {ansi(31, 'π helix')} have symbols {ansi(36, 'G')}, {ansi(36, 'H')} and {ansi(36, 'I')} and are recognized \n"
              f"  by having a repetitive sequence of {ansi(31, 'hydrogen bonds')} in which the residues are three,\n"
              f"  four, or five residues apart respectively. Two types of beta sheet structures exist;\n"
              f"  a beta bridge has symbol {ansi(36, 'B')} while longer sets of {ansi(31, 'hydrogen bonds')} and {ansi(31, 'beta bulges')} have symbol {ansi(36, 'E')}.\n"
              f"  {ansi(36, 'T')} is used for {ansi(31, 'turns')}, featuring hydrogen bonds typical of helices,\n"
              f"  {ansi(36, 'S')} is used for regions of high curvature, and a blank (or space) is used if no other rule applies,\n"
              f"  referring to {ansi(31, 'loops')}. These eight types are usually grouped into three larger classes:\n"
              f"  {ansi(31, 'helix')} ({ansi(36, 'G')}, {ansi(36, 'H')} and {ansi(36, 'I')}), {ansi(31, 'strand')} ({ansi(36, 'E')}"
              f" and {ansi(36, 'B')}) and {ansi(31, 'loop')} ({ansi(36, 'S')}, {ansi(36, 'T')}, and {ansi(36, 'C')}, where {ansi(36, 'C')}"
              f" sometimes is represented also as blank space).\n"
              f"( From Wikipedia [https://en.wikipedia.org/wiki/DSSP_(hydrogen_bond_estimation_algorithm)] )\n\n")

        if usr == "y":
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            Saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if Saving_path == "def":
                Saving_path = f"/home/{nameOfUser}/Documents/Test3"

            print("We are now OPERATING on following path:", os.getcwd())  # The current path
            print(" ")

            source_path = input(f"What PATH do you in mind\
                                (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2"

            os.chdir(source_path)
            source_folder_name = os.popen("ls").read()[:-1]

            print(f"{source_path}/{source_folder_name}")

            os.chdir(f"{source_path}/{source_folder_name}")
            proteins_path = os.popen("ls").readlines()
            for i in range(len(proteins_path)):
                proteins_path[i] = proteins_path[i][:-1]

            number_of_protein = 0
            for protein_in_path in proteins_path:
                ###########################################################################################################
                # The main work of Analysis and gromacs commands.
                # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                print("We are now OPERATING on following path:", os.getcwd())

                print(" ")  # Second folders name list
                print(" ")
                print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                os.popen(f'ls').readlines()

                # Finding Inputs(.tpr)
                files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_tpr)

                print(" ")
                print(" ")
                list_files_root_tpr.sort()
                print(list_files_root_tpr)
                print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                # Finding Inputs(.xtc)
                files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_xtc)

                print(" ")
                print(" ")
                list_files_root_xtc.sort()
                print(list_files_root_xtc)
                print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                print(" ")
                print(" ")
                c = 0
                index = []
                while c <= len(list_files_root_tpr):  # Making folders
                    os.chdir(Saving_path)
                    c += 1
                    if c == len(list_files_root_tpr) + 1:
                        break
                    os.system(f"mkdir do_dssp_{c}")
                    os.chdir(f"{Saving_path}/do_dssp_{c}")
                    index.append(f"do_dssp_{c}")
                print(index)

                print(" ")
                print(" ")
                each_root_tpr = files_roots_tpr.split()  # Copying the inputs to eac folders that we made(.tpr)
                print(each_root_tpr)
                length = len(each_root_tpr)
                i = 0
                list_tpr = []
                while i <= length - 1:
                    list_roots = each_root_tpr[i]  # .tpr
                    print(list_roots[1:])
                    print(" ")
                    print(" ")
                    each_folder = each_root_tpr[i].split("/")
                    list_tpr.append(each_folder[-1])
                    print(each_folder)
                    print(" ")
                    print(" ")
                    print(len(each_folder[-1]))
                    print(" ")
                    print(" ")
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    i += 1

                print(" ")
                print(" ")
                each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                print(each_root_xtc)
                length2 = len(each_root_xtc)
                i = 0
                list_xtc = []

                The_end_of_the_bar = (len(list_files_root_tpr)) * 6
                counting = 0

                while i <= length2 - 1:
                    list_roots = each_root_xtc[i]  # .xtc
                    print(list_roots[1:])
                    each_folder = each_root_xtc[i].split("/")
                    list_xtc.append(each_folder[-1])
                    print(each_folder)
                    print(len(each_folder[-1]))
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    os.chdir(Saving_path + "/" + index[i])

                    gmx = gmxExe
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump0_200-map -o ss0_200.xpm -sc scount0_200.xvg -tu ps -b 0 -e 200")  # Operating
                    os.system(f"{gmx} xpm2ps -f ss0_200.xpm -o ss0_200.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump200_400-map -o ss200_400.xpm -sc scount200_400.xvg -tu ps -b 200 -e 400")
                    os.system(f"{gmx} xpm2ps -f ss200_400.xpm -o ss200_400.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump400_600-map -o ss400_600.xpm -sc scount400_600.xvg -tu ps -b 400 -e 600")
                    os.system(f"{gmx} xpm2ps -f ss400_600.xpm -o ss400_600.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump600_800-map -o ss600_800.xpm -sc scount600_800.xvg -tu ps -b 600 -e 800")
                    os.system(f"{gmx} xpm2ps -f ss600_800.xpm -o ss600_800.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump800_1000-map -o ss800_1000.xpm -sc scount800_1000.xvg -tu ps -b 800 -e 1000")
                    os.system(f"{gmx} xpm2ps -f ss800_1000.xpm -o ss800_1000.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)

                    # The 950 to 1000 ps of simulation
                    os.system(
                        f"echo 1 | {gmx} do_dssp -f {list_xtc[i]} -s {list_tpr[i]} -ssdump ssdump950_1000-map -o ss950_1000.xpm -sc scount950_1000.xvg -tu ps -b 950 -e 1000")
                    os.system(f"{gmx} xpm2ps -f ss950_1000.xpm -o ss950_1000.eps")
                    counting += 1
                    The_Processor(The_end_of_the_bar, counting)

                    i += 1

                print(" ")
                print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                ############################################################################################################
                # Reading the .dat file for measuring the number of each dssp type of structure in each residues and make a excell file of these datas.
                paths = []
                path = Saving_path
                os.chdir(path)
                paths.append(path)
                ans = len(os.popen("ls").readlines())
                print(ans)

                for i in range(1, int(ans) + 1):

                    os.chdir(path + "/" + f"do_dssp_{i}")

                    print(os.getcwd())

                    list_files_root_dssp = os.popen(f'find -name "dssp*"').readlines()
                    print(" ")
                    print(" ")
                    print(list_files_root_dssp)
                    print(f"The Number of files(Inputs): {len(list_files_root_dssp)}")

                    i = 0
                    while i <= len(list_files_root_dssp) - 1:
                        list_files_root_dssp[i] = list_files_root_dssp[i][2:-1]
                        i += 1
                    print(list_files_root_dssp)

                    E_l = []
                    H_l = []
                    B_l = []
                    G_l = []
                    I_l = []
                    T_l = []
                    S_l = []
                    UD_l = []
                    residue_l = []
                    index = []

                    DSSPTypes = ["3-Helix", "5-Helix", "A-Helix", "Turn", "Bend", "B-Bridge", "B-Sheet",
                                 "Coil"]  # The set of DSSPTypes for checking
                    DSSPCodes = ["G", "I", "H", "T", "S", "B", "E", "UD"]

                    for i in ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"]:
                        print(f"----------------------------------{i}----------------------------------")
                        dssp_count = open(f"./scount{i}.xvg").readlines()

                        checkList = {}  # The dictionary for finding the '@ s1' line
                        linecount_checkList = 0  # Start from zero
                        for j in range(len(dssp_count)):
                            checkList[dssp_count[j][:4]] = linecount_checkList
                            linecount_checkList += 1
                        print(checkList)
                        starterLegend = checkList["@ s1"]
                        print(starterLegend)

                        # Finding the number of legends from s1 to sn:
                        True_False = True
                        start = starterLegend
                        legendsDic = {}
                        dsspTypes = []
                        while True_False == True:
                            c = dssp_count[start][:3]
                            print(c)

                            if c != "@ s":
                                True_False = False
                            else:
                                print(dssp_count[start][13:-2])
                                dsspTypes.append(dssp_count[start][13:-2])  # The Dssp types that .xvg file have
                                legendsDic[dssp_count[start][
                                           13:-2]] = []  # Making a list in a dic for each legends that the key is the The number of column that the data is
                                legendsDic[dssp_count[start][13:-2]] \
                                    .append(int(dssp_count[start][
                                                    3]))  # first constituent(with 0 index): The name of each legends
                                legendsDic[dssp_count[start][13:-2]] \
                                    .append([])  # second constituent(with 1 index): a list for the values of data
                                start += 1

                        print(f"The legend dictionary: {legendsDic}")

                        print(f"The end-line of legends: {start}")

                        HList = []
                        BList = []
                        EList = []
                        GList = []
                        UDList = []
                        TList = []
                        SList = []
                        IList = []
                        dsspTypesNotHere = []
                        for j in range(start, (
                                len(dssp_count) - 2)):  # The values for each time (10ps to 10ps) (The second constituent)
                            for typeName in DSSPTypes:
                                exitOrNot = dsspTypes.count(typeName)
                                if exitOrNot == 0:
                                    dsspTypesNotHere.append(typeName)  # The Dssp Types that are not here in .xvg file
                                    pass
                                else:
                                    columnNum = legendsDic[typeName][0]
                                    columnNum -= 1
                                    legendsDic[typeName][1].append(
                                        int(dssp_count[j][14 + (columnNum * 6):14 + ((columnNum + 1) * 6)]))

                        print(legendsDic)
                        print((len(dssp_count) - 2 - start))
                        print(start)
                        print(len(dssp_count) - 2)
                        for j in range(len(DSSPTypes)):
                            exitOrNot = dsspTypes.count(DSSPTypes[j])
                            if exitOrNot == 0:
                                pass
                            else:
                                legendsDic[DSSPTypes[j]] = (sum(legendsDic[DSSPTypes[j]][1])) / (
                                        len(dssp_count) - 2 - start)  # The second constituent become The average

                        if dsspTypesNotHere != []:
                            for dt in range(len(dsspTypesNotHere)):  # Adding the Not found Dssp types to the legendsDic
                                legendsDic[dsspTypesNotHere[dt]] = 0

                        for k in range(
                                len(DSSPTypes)):  # Changing the name of the keys to DSSP codes(H, B, E G, T, S, UD)
                            value = legendsDic.pop(DSSPTypes[k])
                            legendsDic[DSSPCodes[k]] = value

                        legendsDic["Res"] = round(
                            legendsDic["H"] + legendsDic["B"] + legendsDic["G"] + legendsDic["E"] + legendsDic["T"] +
                            legendsDic["S"] + legendsDic["UD"] + legendsDic["I"])

                        print(legendsDic)

                        E_l.append(legendsDic["E"])
                        H_l.append(legendsDic["H"])
                        B_l.append(legendsDic["B"])
                        G_l.append(legendsDic["G"])
                        T_l.append(legendsDic["T"])
                        S_l.append(legendsDic["S"])
                        UD_l.append(legendsDic["UD"])
                        I_l.append(legendsDic["I"])
                        residue_l.append(legendsDic["Res"])

                        index.append(i + "(ns)")

                    dic_dssp = {"H": H_l, "B": B_l, \
                                "E": E_l, "G": G_l, "I": I_l, "T": T_l, "S": S_l, "UD": UD_l, "Res": residue_l}
                    print("DSSP:")
                    print(dic_dssp)
                    print(index)

                    n = pd.DataFrame(data=dic_dssp, index=index)
                    print(n)

                    writer = pd.ExcelWriter('The_DSSP_analysis.xlsx', engine='xlsxwriter')
                    n.to_excel(writer, sheet_name="The DSSP Analysis")
                    writer.save()

                ############################################################################################################
                # Changing the name of the folders of each Acceleration simulation's analysis
                ans = Saving_path
                dir = ans
                os.chdir(dir)
                folders = os.popen("ls").readlines()

                print(folders)

                for i in folders:
                    os.chdir(f"{dir}/do_dssp_{i[8:-1]}")
                    file = os.popen("ls").readlines()
                    print(file)
                    acc_num = file[6][7:-5]
                    print(acc_num + "<")
                    if acc_num == "":
                        acc_num = 0
                    os.chdir(f"{dir}")
                    os.system(f"mv do_dssp_{i[8:-1]} do_dssp_{acc_num}")

                ############################################################################################################
                # Making a excel folder for the whole Analysis of each Acceleration simulation and partitioning it.
                ans2 = Saving_path
                os.chdir(ans2)

                lis_folders = os.popen("ls").readlines()
                print(lis_folders)

                H0_1000 = []
                B0_1000 = []
                E0_1000 = []
                G0_1000 = []
                I0_1000 = []
                T0_1000 = []
                S0_1000 = []
                UD0_1000 = []
                Total_Helix0_1000 = []
                Total_Sheet0_1000 = []
                Total_Loop0_1000 = []
                Total_SS0_1000 = []
                Dic_DSSP0_1000 = {}

                H800_1000 = []
                B800_1000 = []
                E800_1000 = []
                G800_1000 = []
                I800_1000 = []
                T800_1000 = []
                S800_1000 = []
                UD800_1000 = []
                Total_Helix800_1000 = []
                Total_Sheet800_1000 = []
                Total_Loop800_1000 = []
                Total_SS800_1000 = []
                Dic_DSSP800_1000 = {}

                H950_1000 = []
                B950_1000 = []
                E950_1000 = []
                G950_1000 = []
                I950_1000 = []
                T950_1000 = []
                S950_1000 = []
                UD950_1000 = []
                Total_Helix950_1000 = []
                Total_Sheet950_1000 = []
                Total_Loop950_1000 = []
                Total_SS950_1000 = []
                Dic_DSSP950_1000 = {}

                index = []

                for i in range(len(lis_folders)):
                    os.chdir(f"{ans2}/{lis_folders[i][:-1]}")
                    lis_files = os.popen("ls").readlines()
                    for j in range(len(lis_files)):
                        lis_files[j] = lis_files[j][:-1]
                    print(lis_files)

                    datas = pd.read_excel("The_DSSP_analysis.xlsx")

                    c0_1000 = 0
                    sumH0_1000 = 0
                    sumB0_1000 = 0
                    sumE0_1000 = 0
                    sumG0_1000 = 0
                    sumI0_1000 = 0
                    sumT0_1000 = 0
                    sumS0_1000 = 0
                    sumUD0_1000 = 0
                    sumTotal_Helix0_1000 = 0
                    sumTotal_Sheet0_1000 = 0
                    sumTotal_Loop0_1000 = 0
                    sumTotal_SS0_1000 = 0
                    while c0_1000 != 5:
                        sumH0_1000 += datas["H"][c0_1000]
                        sumB0_1000 += datas["B"][c0_1000]
                        sumE0_1000 += datas["E"][c0_1000]
                        sumG0_1000 += datas["G"][c0_1000]
                        sumI0_1000 += datas["I"][c0_1000]
                        sumT0_1000 += datas["T"][c0_1000]
                        sumS0_1000 += datas["S"][c0_1000]
                        sumUD0_1000 += datas["UD"][c0_1000]
                        sumTotal_Helix0_1000 += (datas["H"][c0_1000] + datas["G"][c0_1000] + datas["I"][c0_1000])
                        sumTotal_Sheet0_1000 += (datas["B"][c0_1000] + datas["E"][c0_1000])
                        sumTotal_Loop0_1000 += (datas["T"][c0_1000] + datas["S"][c0_1000] + datas["UD"][c0_1000])
                        sumTotal_SS0_1000 += (datas["H"][c0_1000] +
                                              datas["B"][c0_1000] +
                                              datas["E"][c0_1000] +
                                              datas["G"][c0_1000] +
                                              datas["I"][c0_1000])
                        c0_1000 += 1
                    aveH0_1000 = sumH0_1000 / 100
                    H0_1000.append(aveH0_1000)
                    aveB0_1000 = sumB0_1000 / 100
                    B0_1000.append(aveB0_1000)
                    aveE0_1000 = sumE0_1000 / 100
                    E0_1000.append(aveE0_1000)
                    aveG0_1000 = sumG0_1000 / 100
                    G0_1000.append(aveG0_1000)
                    aveI0_1000 = sumI0_1000 / 100
                    I0_1000.append(aveI0_1000)
                    aveT0_1000 = sumT0_1000 / 100
                    T0_1000.append(aveT0_1000)
                    aveS0_1000 = sumS0_1000 / 100
                    S0_1000.append(aveS0_1000)
                    aveUD0_1000 = sumUD0_1000 / 100
                    UD0_1000.append(aveUD0_1000)
                    aveTotal_Helix0_1000 = sumTotal_Helix0_1000 / 100
                    Total_Helix0_1000.append(aveTotal_Helix0_1000)
                    aveTotal_Sheet0_1000 = sumTotal_Sheet0_1000 / 100
                    Total_Sheet0_1000.append(aveTotal_Sheet0_1000)
                    aveTotal_Loop0_1000 = sumTotal_Loop0_1000 / 100
                    Total_Loop0_1000.append(aveTotal_Loop0_1000)
                    aveTotal_SS0_1000 = sumTotal_SS0_1000 / 100
                    Total_SS0_1000.append(aveTotal_SS0_1000)

                    aveH800_1000 = datas["H"][4]
                    H800_1000.append(aveH800_1000)
                    aveB800_1000 = datas["B"][4]
                    B800_1000.append(aveB800_1000)
                    aveE800_1000 = datas["E"][4]
                    E800_1000.append(aveE800_1000)
                    aveG800_1000 = datas["G"][4]
                    G800_1000.append(aveG800_1000)
                    aveI800_1000 = datas["I"][4]
                    I800_1000.append(aveI800_1000)
                    aveT800_1000 = datas["T"][4]
                    T800_1000.append(aveT800_1000)
                    aveS800_1000 = datas["S"][4]
                    S800_1000.append(aveS800_1000)
                    aveUD800_1000 = datas["UD"][4]
                    UD800_1000.append(aveUD800_1000)
                    aveTotal_Helix800_1000 = (datas["H"][4] + datas["G"][4] + datas["I"][4])
                    Total_Helix800_1000.append(aveTotal_Helix800_1000)
                    aveTotal_Sheet800_1000 = (datas["B"][4] + datas["E"][4])
                    Total_Sheet800_1000.append(aveTotal_Sheet800_1000)
                    aveTotal_Loop800_1000 = (datas["S"][4] + datas["T"][4] + datas["UD"][4])
                    Total_Loop800_1000.append(aveTotal_Loop800_1000)
                    aveTotal_SS800_1000 = (datas["H"][4] +
                                           datas["B"][4] +
                                           datas["E"][4] +
                                           datas["G"][4] +
                                           datas["I"][4])
                    Total_SS800_1000.append(aveTotal_SS800_1000)

                    aveH950_1000 = datas["H"][5]
                    H950_1000.append(aveH950_1000)
                    aveB950_1000 = datas["B"][5]
                    B950_1000.append(aveB950_1000)
                    aveE950_1000 = datas["E"][5]
                    E950_1000.append(aveE950_1000)
                    aveG950_1000 = datas["G"][5]
                    G950_1000.append(aveG950_1000)
                    aveI950_1000 = datas["I"][5]
                    I950_1000.append(aveI950_1000)
                    aveT950_1000 = datas["T"][5]
                    T950_1000.append(aveT950_1000)
                    aveS950_1000 = datas["S"][5]
                    S950_1000.append(aveS950_1000)
                    aveUD950_1000 = datas["UD"][5]
                    UD950_1000.append(aveUD950_1000)
                    aveTotal_Helix950_1000 = (datas["H"][5] + datas["G"][5] + datas["I"][5])
                    Total_Helix950_1000.append(aveTotal_Helix950_1000)
                    aveTotal_Sheet950_1000 = (datas["B"][5] + datas["E"][5])
                    Total_Sheet950_1000.append(aveTotal_Sheet950_1000)
                    aveTotal_Loop950_1000 = (datas["S"][5] + datas["T"][5] + datas["UD"][5])
                    Total_Loop950_1000.append(aveTotal_Loop950_1000)
                    aveTotal_SS950_1000 = (datas["H"][5] +
                                           datas["B"][5] +
                                           datas["E"][5] +
                                           datas["G"][5] +
                                           datas["I"][5])
                    Total_SS950_1000.append(aveTotal_SS950_1000)

                    index.append(lis_folders[i][8:-1])

                    Dic_DSSP0_1000["α-helix"] = H0_1000
                    Dic_DSSP0_1000["residue in isolated β-bridge"] = B0_1000
                    Dic_DSSP0_1000["extended strand, participates in β ladder"] = E0_1000
                    Dic_DSSP0_1000["3-helix (310 helix)"] = G0_1000
                    Dic_DSSP0_1000["5 helix (π-helix)"] = I0_1000
                    Dic_DSSP0_1000["hydrogen bonded turn"] = T0_1000
                    Dic_DSSP0_1000["bend"] = S0_1000
                    Dic_DSSP0_1000["random coil"] = UD0_1000
                    Dic_DSSP0_1000["total helices"] = Total_Helix0_1000
                    Dic_DSSP0_1000["total strands"] = Total_Sheet0_1000
                    Dic_DSSP0_1000["total loops"] = Total_Loop0_1000
                    Dic_DSSP0_1000["total regular structure"] = Total_SS0_1000

                    Dic_DSSP800_1000["α-helix"] = H800_1000
                    Dic_DSSP800_1000["residue in isolated β-bridge"] = B800_1000
                    Dic_DSSP800_1000["extended strand, participates in β ladder"] = E800_1000
                    Dic_DSSP800_1000["3-helix (310 helix)"] = G800_1000
                    Dic_DSSP800_1000["5 helix (π-helix)"] = I800_1000
                    Dic_DSSP800_1000["hydrogen bonded turn"] = T800_1000
                    Dic_DSSP800_1000["bend"] = S800_1000
                    Dic_DSSP800_1000["random coil"] = UD800_1000
                    Dic_DSSP800_1000["total helices"] = Total_Helix800_1000
                    Dic_DSSP800_1000["total strands"] = Total_Sheet800_1000
                    Dic_DSSP800_1000["total loops"] = Total_Loop800_1000
                    Dic_DSSP800_1000["total regular structure"] = Total_SS800_1000

                    Dic_DSSP950_1000["α-helix"] = H950_1000
                    Dic_DSSP950_1000["residue in isolated β-bridge"] = B950_1000
                    Dic_DSSP950_1000["extended strand, participates in β ladder"] = E950_1000
                    Dic_DSSP950_1000["3-helix (310 helix)"] = G950_1000
                    Dic_DSSP950_1000["5 helix (π-helix)"] = I950_1000
                    Dic_DSSP950_1000["hydrogen bonded turn"] = T950_1000
                    Dic_DSSP950_1000["bend"] = S950_1000
                    Dic_DSSP950_1000["random coil"] = UD950_1000
                    Dic_DSSP950_1000["total helices"] = Total_Helix950_1000
                    Dic_DSSP950_1000["total strands"] = Total_Sheet950_1000
                    Dic_DSSP950_1000["total loops"] = Total_Loop950_1000
                    Dic_DSSP950_1000["total regular structure"] = Total_SS950_1000

                    print(Dic_DSSP0_1000)
                    print(Dic_DSSP800_1000)
                    print(Dic_DSSP950_1000)

                os.chdir(ans2)

                dic0_1000 = pd.DataFrame(data=Dic_DSSP0_1000, index=index)
                print(dic0_1000)
                writer1 = pd.ExcelWriter('Total_DSSP.xlsx', engine='xlsxwriter')
                dic0_1000.to_excel(writer1, sheet_name="The DSSP Analysis")
                writer1.save()

                dic800_1000 = pd.DataFrame(data=Dic_DSSP800_1000, index=index)
                writer2 = pd.ExcelWriter('800_1000_DSSP.xlsx', engine='xlsxwriter')
                dic800_1000.to_excel(writer2, sheet_name="The DSSP Analysis")
                writer2.save()

                dic950_1000 = pd.DataFrame(data=Dic_DSSP950_1000, index=index)
                print("Hello:\n", dic950_1000)
                writer3 = pd.ExcelWriter('950_1000_DSSP.xlsx', engine='xlsxwriter')
                dic950_1000.to_excel(writer3, sheet_name="The DSSP Analysis")
                writer3.save()
                ############################################################################################################
                # Making a folder to moving the files and folders in that folder.
                os.chdir(f"{source_path}/{source_folder_name}")
                pro_name = os.popen("ls").readlines()
                pro_name = pro_name[number_of_protein]
                for i in range(len(pro_name)):
                    if pro_name[i] == "_":
                        pro_name1 = pro_name[0:i]
                        pro_name2 = pro_name[i + 1:-1]
                        break

                os.chdir(Saving_path)
                files = os.popen("ls").readlines()
                for i in range(len(files)):
                    files[i] = files[i][:-1]

                Analysis_folders = pro_name1 + "_DSSP_" + pro_name2
                os.system(f"mkdir {Analysis_folders}")

                for i in files:
                    print(f"Hello{i}")
                    print(f"{Analysis_folders}")
                    os.system(f"mv {i} {Analysis_folders}")

                print(f"={Analysis_folders}")
                print(type(Analysis_folders))
                print(Analysis_folders == "1fsv_DSSP_0to1ns_2")
                print(f"={final_Saving_path}")
                print(type(final_Saving_path))
                print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                print(f"mv {Analysis_folders} {final_Saving_path}")
                os.system(f"mv {Analysis_folders} {final_Saving_path}")

                number_of_protein += 1
                ############################################################################################################

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style=3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style=7)} to exit the program) \n")
    # -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
    elif alys_name == "rmsd":                                                                                           # The RMSD analysis
        ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'PLEASE READ', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        if usr == "y":
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            Saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if Saving_path == "def":
                Saving_path = f"/home/{nameOfUser}/Documents/Test3"

            print("We are now OPERATING on following path:", os.getcwd())  # The current path
            print(" ")

            source_path = input(f"What PATH do you in mind\
                                (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2"

            os.chdir(source_path)
            source_folder_name = os.popen("ls").read()[:-1]

            os.chdir(f"{source_path}/{source_folder_name}")
            proteins_path = os.popen("ls").readlines()
            for i in range(len(proteins_path)):
                proteins_path[i] = proteins_path[i][:-1]

            number_of_protein = 0
            for protein_in_path in proteins_path:
                ###########################################################################################################
                # The main work of Analysis and gromacs commands.
                # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                print("We are now OPERATING on following path:", os.getcwd())

                print(" ")  # Second folders name list
                print(" ")
                print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                os.popen(f'ls').readlines()

                # Finding Inputs(.tpr)
                files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_tpr)

                print(" ")
                print(" ")
                list_files_root_tpr.sort()
                print(list_files_root_tpr)
                print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                # Finding Inputs(.xtc)
                files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_xtc)

                print(" ")
                print(" ")
                list_files_root_xtc.sort()
                print(list_files_root_xtc)
                print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                print(" ")
                print(" ")
                c = 0
                index = []
                while c <= len(list_files_root_tpr):  # Making folders
                    os.chdir(Saving_path)
                    c += 1
                    if c == len(list_files_root_tpr) + 1:
                        break
                    os.system(f"mkdir rms_{c}")
                    os.chdir(f"{Saving_path}/rms_{c}")
                    index.append(f"rms_{c}")
                print(index)

                print(" ")
                print(" ")
                each_root_tpr = files_roots_tpr.split()  # Copying the inputs to each folders that we made(.tpr)
                print(each_root_tpr)
                length = len(each_root_tpr)
                i = 0
                list_tpr = []
                while i <= length - 1:
                    list_roots = each_root_tpr[i]  # .tpr
                    print(list_roots[1:])
                    print(" ")
                    print(" ")
                    each_folder = each_root_tpr[i].split("/")
                    list_tpr.append(each_folder[-1])
                    print(each_folder)
                    print(" ")
                    print(" ")
                    print(len(each_folder[-1]))
                    print(" ")
                    print(" ")
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    i += 1

                print(" ")
                print(" ")
                each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                print(each_root_xtc)
                length2 = len(each_root_xtc)
                i = 0
                list_xtc = []

                The_end_of_the_bar = (len(list_files_root_tpr)) * 5
                counting = 0

                while i <= length2 - 1:
                    list_roots = each_root_xtc[i]  # .xtc
                    print(list_roots[1:])
                    each_folder = each_root_xtc[i].split("/")
                    list_xtc.append(each_folder[-1])
                    print(each_folder)
                    print(len(each_folder[-1]))
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    os.chdir(Saving_path + "/" + index[i])

                    gmx = gmxExe
                    os.system(
                        f"echo 1 1 | {gmx} rms -f {list_xtc[i]} -s {list_tpr[i]} -o rms.xvg")  # Operating

                    counting += 1

                    The_Processor(len(list_files_root_tpr), counting)

                    rms = open(f"rms.xvg").readlines()
                    length12 = len(rms)

                    time_list = []
                    rms_list = []
                    index_1 = []
                    c = 0

                    for k in range(18, length12):

                        time_list.append((float(rms[k][0:12])))

                        rms_list.append(float(rms[k][16:25]))

                        c += 1
                        index_1.append(c)

                    ex_dict = {"time": time_list, "RMSD": rms_list}

                    n = pd.DataFrame(data=ex_dict, index=index_1)

                    writer = pd.ExcelWriter('The_analysis_RMSD.xlsx', engine='xlsxwriter')
                    n.to_excel(writer, sheet_name="The RMSD Analysis")
                    writer.save()

                    i += 1
                    print(" ")
                    print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                ############################################################################################################
                # Changing the name of the folders of each Acceleration simulation's analysis
                ans = Saving_path
                dir = ans
                os.chdir(dir)
                folders = os.popen("ls").readlines()

                print(folders)

                for i in folders:
                    os.chdir(f"{dir}/rms{i[3:-1]}")
                    file = os.popen("ls").readlines()
                    print(file)
                    acc_num = file[0][7:-5]  # the first number is the number of files in each folders.
                    print(acc_num + "<")
                    if acc_num == "":
                        acc_num = 0
                    os.chdir(f"{dir}")
                    os.system(f"mv rms{i[3:-1]} rms_{acc_num}")

                ############################################################################################################
                # Making a excel folder for the whole Analysis of each Acceleration simulation and partitioning it.
                ans2 = Saving_path
                os.chdir(ans2)

                lis_folders = os.popen("ls").readlines()
                print(lis_folders)
                dic_datas0_1000 = {}
                dic_datas800_1000 = {}
                dic_data_AccToTime = {}
                for i in range(len(lis_folders)):

                    os.chdir(f"{ans2}/{lis_folders[i][:-1]}")
                    lis_files = os.popen("ls").readlines()
                    for j in range(len(lis_files)):
                        lis_files[j] = lis_files[j][:-1]
                    print(lis_files)
                    name = lis_files[1][-8:-4]

                    datas = pd.read_excel("The_analysis_RMSD.xlsx")
                    c0_1000 = 0
                    sum0_1000 = 0
                    while c0_1000 != len(datas["RMSD"]):
                        sum0_1000 += datas["RMSD"][c0_1000]
                        c0_1000 += 1
                    ave0_1000 = sum0_1000 / 100

                    dic_datas0_1000[lis_folders[i][4:-1]] = ave0_1000

                    c800_1000 = 80
                    sum800_1000 = 0
                    while c800_1000 != len(datas["RMSD"]):
                        sum800_1000 += datas["RMSD"][c800_1000]
                        c800_1000 += 1
                    ave800_1000 = sum800_1000 / 20

                    dic_datas800_1000[lis_folders[i][4:-1]] = ave800_1000

                    dic_data_AccToTime[f"{name}_RMSD"] = datas["RMSD"]

                print(dic_datas0_1000)
                print(dic_datas800_1000)
                print(dic_data_AccToTime)

                dic_data_AccToTime["time"] = datas["time"]

                os.chdir(ans2)

                dic0_1000 = pd.DataFrame(data=dic_datas0_1000, index=[1])
                print(dic0_1000)
                writer1 = pd.ExcelWriter('Total_RMSD.xlsx', engine='xlsxwriter')
                dic0_1000.to_excel(writer1, sheet_name="The RMSD Analysis")
                writer1.save()

                dic800_1000 = pd.DataFrame(data=dic_datas800_1000, index=[1])
                writer2 = pd.ExcelWriter('800_1000_RMSD.xlsx', engine='xlsxwriter')
                dic800_1000.to_excel(writer2, sheet_name="The RMSD Analysis")
                writer2.save()

                dicAccToTime = pd.DataFrame(data=dic_data_AccToTime)
                writer3 = pd.ExcelWriter('Acc_Time_RMSD.xlsx', engine='xlsxwriter')
                dicAccToTime.to_excel(writer3, sheet_name="The RMSD Analysis")
                writer3.save()
                ############################################################################################################
                # Making a folder to moving the files and folders in that folder.
                os.chdir(f"{source_path}/{source_folder_name}")
                pro_name = os.popen("ls").readlines()
                pro_name = pro_name[number_of_protein]
                for i in range(len(pro_name)):
                    if pro_name[i] == "_":
                        pro_name1 = pro_name[0:i]
                        pro_name2 = pro_name[i + 1:-1]
                        break

                os.chdir(Saving_path)
                files = os.popen("ls").readlines()
                for i in range(len(files)):
                    files[i] = files[i][:-1]

                Analysis_folders = pro_name1 + "_RMSD_" + pro_name2
                os.system(f"mkdir {Analysis_folders}")

                for i in files:
                    print(f"Hello{i}")
                    print(f"{Analysis_folders}")
                    os.system(f"mv {i} {Analysis_folders}")

                print(f"={Analysis_folders}")
                print(type(Analysis_folders))
                print(Analysis_folders == "1fsv_RMSD_0to1ns_2")
                print(f"={final_Saving_path}")
                print(type(final_Saving_path))
                print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                print(f"mv {Analysis_folders} {final_Saving_path}")
                os.system(f"mv {Analysis_folders} {final_Saving_path}")

                number_of_protein += 1
                ############################################################################################################

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style = 3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style = 7)} to exit the program) \n")

    elif alys_name == "rmsf":                                                                                         # The RMSF analysis
                ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'PLEASE READ', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        if usr == "y":
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            Saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if Saving_path == "def":
                Saving_path = f"/home/{nameOfUser}/Documents/Test3"

            print("We are now OPERATING on following path:", os.getcwd())  # The current path
            print(" ")

            source_path = input(f"What PATH do you in mind\
                                (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2"

            # Choosing the "rmsf_atomOrRes" or rmsf_10frame:
            ans3 = input(f"What do you want to {ansi(32, 'do', style=3)} with {ansi(33, 'RMSF', style=0)}?\n(choose between '{ansi(94, '10', style=3)}' for 10 frames and '{ansi(35, 'ar', style=3)}' for atom or residue rmsf)")
            if ans3 == "10":
                rmsf_type = "10_frames"
            elif ans3 == "ar":
                rmsf_type = "atomOrResidue"
                ans4 = input("Do you wish doing RMSD over Atoms or Residues?(a/r)\n")

            os.chdir(source_path)
            source_folder_name = os.popen("ls").read()[:-1]

            os.chdir(f"{source_path}/{source_folder_name}")
            proteins_path = os.popen("ls").readlines()
            for i in range(len(proteins_path)):
                proteins_path[i] = proteins_path[i][:-1]

            number_of_protein = 0
            for protein_in_path in proteins_path:
                if ans3 == "ar":

                    ###########################################################################################################
                    # The main work of Analysis and gromacs commands.
                    # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                    os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                    print("We are now OPERATING on following path:", os.getcwd())

                    print(" ")  # Second folders name list
                    print(" ")
                    print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                    os.popen(f'ls').readlines()

                    # Finding Inputs(.tpr)
                    files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                    list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                    print(" ")
                    print(" ")
                    print("The ROOT of your wanted files: ")
                    print(files_roots_tpr)

                    print(" ")
                    print(" ")
                    list_files_root_tpr.sort()
                    print(list_files_root_tpr)
                    print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                    # Finding Inputs(.xtc)
                    files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                    list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                    print(" ")
                    print(" ")
                    print("The ROOT of your wanted files: ")
                    print(files_roots_xtc)

                    print(" ")
                    print(" ")
                    list_files_root_xtc.sort()
                    print(list_files_root_xtc)
                    print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                    print(" ")
                    print(" ")
                    c = 0
                    index = []
                    while c <= len(list_files_root_tpr):  # Making folders
                        os.chdir(Saving_path)
                        c += 1
                        if c == len(list_files_root_tpr) + 1:
                            break
                        os.system(f"mkdir rmsf_{rmsf_type}_{c}")
                        os.chdir(f"{Saving_path}/rmsf_{rmsf_type}_{c}")
                        index.append(f"rmsf_{rmsf_type}_{c}")
                    print(index)

                    print(" ")
                    print(" ")
                    each_root_tpr = files_roots_tpr.split()  # Copying the inputs to eac folders that we made(.tpr)
                    print(each_root_tpr)
                    length = len(each_root_tpr)
                    i = 0
                    list_tpr = []
                    while i <= length - 1:
                        list_roots = each_root_tpr[i]  # .tpr
                        print(list_roots[1:])
                        print(" ")
                        print(" ")
                        each_folder = each_root_tpr[i].split("/")
                        list_tpr.append(each_folder[-1])
                        print(each_folder)
                        print(" ")
                        print(" ")
                        print(len(each_folder[-1]))
                        print(" ")
                        print(" ")
                        os.chdir(
                            f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                        os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                        print(Saving_path + "/" + index[i])

                        i += 1

                    print(" ")
                    print(" ")
                    each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                    print(each_root_xtc)
                    length2 = len(each_root_xtc)
                    i = 0
                    list_xtc = []

                    The_end_of_the_bar = (len(list_files_root_tpr)) * 5
                    counting = 0

                    while i <= length2 - 1:
                        list_roots = each_root_xtc[i]  # .xtc
                        print(list_roots[1:])
                        each_folder = each_root_xtc[i].split("/")
                        list_xtc.append(each_folder[-1])
                        print(each_folder)
                        print(len(each_folder[-1]))
                        os.chdir(
                            f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                        os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                        print(Saving_path + "/" + index[i])

                        os.chdir(Saving_path + "/" + index[i])

                        gmx = gmxExe

                        if ans4 == "a":
                            os.system(
                                f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf.xvg -ox average.pdb -b 0 -e 1000")  # Operating
                        elif ans4 == "r":
                            os.system(
                                f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf.xvg -ox average.pdb -res -b 0 -e 1000")  # Operating

                        counting += 1

                        The_Processor(len(list_files_root_tpr), counting)
                        os.getcwd()

                        rmsf = open(f"rmsf.xvg").readlines()
                        length12 = len(rmsf)

                        cons_list = []
                        rmsf_list = []
                        index_1 = []
                        c = 0

                        for k in range(17, length12):
                            cons_list.append(int(rmsf[k][:6]))  # residue/atom
                            rmsf_list.append(float(rmsf[k][8:]))  # rmsf
                            c += 1
                            index_1.append(c)

                        ex_dict = {f"{ans4}": cons_list, "rmsf": rmsf_list}

                        n = pd.DataFrame(data=ex_dict, index=index_1)

                        writer = pd.ExcelWriter('The_analysis_rmsf.xlsx', engine='xlsxwriter')
                        n.to_excel(writer, sheet_name="The RMSF Analysis")
                        writer.save()

                        i += 1
                        print(" ")
                        print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                    ############################################################################################################
                    # Changing the name of the folders of each Acceleration simulation's analysis
                    ans = Saving_path
                    dir = ans
                    os.chdir(dir)
                    folders = os.popen("ls").readlines()

                    print(folders)

                    for i in folders:
                        os.chdir(f"{dir}/rmsf_{i[5:-1]}")
                        file = os.popen("ls").readlines()
                        print(file)
                        acc_num = file[1][7:-5]  # the first number is the number of files in each folders.
                        print(acc_num + "<")
                        if acc_num == "":
                            acc_num = 0
                        os.chdir(f"{dir}")
                        os.system(f"mv rmsf_{i[5:-1]} rmsf_{acc_num}")

                    ############################################################################################################
                    # Making a excel folder for the whole Analysis of each Acceleration simulation and partitioning it.
                    ans2 = Saving_path
                    os.chdir(ans2)

                    lis_folders = os.popen("ls").readlines()
                    print(lis_folders)
                    RMSF0_1000 = []
                    dic_datas0_1000 = {}
                    index = []
                    for i in range(len(lis_folders)):

                        os.chdir(f"{ans2}/{lis_folders[i][:-1]}")
                        lis_files = os.popen("ls").readlines()
                        for j in range(len(lis_files)):
                            lis_files[j] = lis_files[j][:-1]
                        print(lis_files)

                        datas = pd.read_excel("The_analysis_rmsf.xlsx")
                        c0_1000 = 0
                        sum0_1000 = 0
                        while c0_1000 != len(cons_list):
                            sum0_1000 += datas["rmsf"][c0_1000]
                            c0_1000 += 1
                        ave0_1000H = sum0_1000 / len(cons_list)

                        RMSF0_1000.append(ave0_1000H)

                        index.append(lis_folders[i][6:-1])

                    dic_datas0_1000["rmsf"] = RMSF0_1000
                    print(dic_datas0_1000)

                    os.chdir(ans2)

                    dic0_1000 = pd.DataFrame(data=dic_datas0_1000, index=index)
                    print(dic0_1000)
                    writer1 = pd.ExcelWriter('Total_rmsf.xlsx', engine='xlsxwriter')
                    dic0_1000.to_excel(writer1, sheet_name="The rmsf Analysis")
                    writer1.save()

                    ############################################################################################################
                    # Making a folder to moving the files and folders in that folder.
                    os.chdir(f"{source_path}/{source_folder_name}")
                    pro_name = os.popen("ls").readlines()
                    pro_name = pro_name[number_of_protein]
                    for i in range(len(pro_name)):
                        if pro_name[i] == "_":
                            pro_name1 = pro_name[0:i]
                            pro_name2 = pro_name[i + 1:-1]
                            break

                    os.chdir(Saving_path)
                    files = os.popen("ls").readlines()
                    for i in range(len(files)):
                        files[i] = files[i][:-1]

                    Analysis_folders = pro_name1 + "_RMSF_" + pro_name2
                    os.system(f"mkdir {Analysis_folders}")

                    for i in files:
                        print(f"Hello{i}")
                        print(f"{Analysis_folders}")
                        os.system(f"mv {i} {Analysis_folders}")

                    print(f"={Analysis_folders}")
                    print(type(Analysis_folders))
                    print(Analysis_folders == "1fsv_RMSF_0to1ns_2")
                    print(f"={final_Saving_path}")
                    print(type(final_Saving_path))
                    print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                    # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                    print(f"mv {Analysis_folders} {final_Saving_path}")
                    os.system(f"mv {Analysis_folders} {final_Saving_path}")

                    number_of_protein += 1

                elif ans3 == "10":
                    ###########################################################################################################
                    # The main work of Analysis and gromacs commands.
                    # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                    os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                    print("We are now OPERATING on following path:", os.getcwd())

                    print(" ")  # Second folders name list
                    print(" ")
                    print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                    os.popen(f'ls').readlines()

                    # Finding Inputs(.tpr)
                    files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                    list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                    print(" ")
                    print(" ")
                    print("The ROOT of your wanted files: ")
                    print(files_roots_tpr)

                    print(" ")
                    print(" ")
                    list_files_root_tpr.sort()
                    print(list_files_root_tpr)
                    print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                    # Finding Inputs(.xtc)
                    files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                    list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                    print(" ")
                    print(" ")
                    print("The ROOT of your wanted files: ")
                    print(files_roots_xtc)

                    print(" ")
                    print(" ")
                    list_files_root_xtc.sort()
                    print(list_files_root_xtc)
                    print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                    print(" ")
                    print(" ")
                    c = 0
                    index = []
                    while c <= len(list_files_root_tpr):  # Making folders
                        os.chdir(Saving_path)
                        c += 1
                        if c == len(list_files_root_tpr) + 1:
                            break
                        os.system(f"mkdir rmsf_{rmsf_type}_{c}")
                        os.chdir(f"{Saving_path}/rmsf_{rmsf_type}_{c}")
                        index.append(f"rmsf_{rmsf_type}_{c}")
                    print(index)

                    print(" ")
                    print(" ")
                    each_root_tpr = files_roots_tpr.split()  # Copying the inputs to eac folders that we made(.tpr)
                    print(each_root_tpr)
                    length = len(each_root_tpr)
                    i = 0
                    list_tpr = []
                    while i <= length - 1:
                        list_roots = each_root_tpr[i]  # .tpr
                        print(list_roots[1:])
                        print(" ")
                        print(" ")
                        each_folder = each_root_tpr[i].split("/")
                        list_tpr.append(each_folder[-1])
                        print(each_folder)
                        print(" ")
                        print(" ")
                        print(len(each_folder[-1]))
                        print(" ")
                        print(" ")
                        os.chdir(
                            f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                        os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                        print(Saving_path + "/" + index[i])

                        i += 1

                    print(" ")
                    print(" ")
                    each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                    print(each_root_xtc)
                    length2 = len(each_root_xtc)
                    i = 0
                    list_xtc = []

                    The_end_of_the_bar = (len(list_files_root_tpr)) * 5
                    counting = 0

                    while i <= length2 - 1:
                        list_roots = each_root_xtc[i]  # .xtc
                        print(list_roots[1:])
                        each_folder = each_root_xtc[i].split("/")
                        list_xtc.append(each_folder[-1])
                        print(each_folder)
                        print(len(each_folder[-1]))
                        os.chdir(
                            f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                        os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                        print(Saving_path + "/" + index[i])

                        os.chdir(Saving_path + "/" + index[i])

                        gmx = gmxExe
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf0_100.xvg -ox average0_100.pdb -b 0 -e 100")  # Operating
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf100_200.xvg -ox average100_200.pdb -b 100 -e 200")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf200_300.xvg -ox average200_300.pdb -b 200 -e 300")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf300_400.xvg -ox average300_400.pdb -b 300 -e 400")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf400_500.xvg -ox average400_500.pdb -b 400 -e 500")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf500_600.xvg -ox average500_600.pdb -b 500 -e 600")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf600_700.xvg -ox average600_700.pdb -b 600 -e 700")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf700_800.xvg -ox average700_800.pdb -b 700 -e 800")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf800_900.xvg -ox average800_900.pdb -b 800 -e 900")
                        os.system(
                            f"echo 1 | {gmx} rmsf -f {list_xtc[i]} -s {list_tpr[i]} -o rmsf900_1000.xvg -ox average900_1000.pdb -b 900 -e 1000")
                        counting += 1

                        The_Processor(len(list_files_root_tpr), counting)
                        os.getcwd()

                        i += 1
                        print(" ")
                        print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                    ############################################################################################################
                    # Changing the name of the folders of each Acceleration simulation's analysis
                    dir = Saving_path
                    os.chdir(dir)
                    folders = os.popen("ls").readlines()

                    print(folders)

                    for i in folders:
                        os.chdir(f"{dir}/rmsf_{i[5:-1]}")
                        file = os.popen("ls").readlines()
                        print(file)
                        acc_num = file[10][7:-5]  # the first number is the number of files in each folders.
                        print(acc_num + "<")
                        if acc_num == "":
                            acc_num = 0
                        os.chdir(f"{dir}")
                        os.system(f"mv rmsf_{i[5:-1]} rmsf_{acc_num}")

                    ############################################################################################################
                    # Pymoling the .pdb files to construct their Pictures.
                    os.chdir(f"{source_path}/{source_folder_name}")     # Name
                    pro_name = os.popen("ls").readlines()
                    pro_name = pro_name[number_of_protein]
                    for i in range(len(pro_name)):
                        if pro_name[i] == "_":
                            pro_name1 = pro_name[0:i]
                            pro_name2 = pro_name[i + 1:-1]
                            break
                    os.chdir(f"{Saving_path}")
                    os.system(f"mkdir .Pics_{pro_name1}_{pro_name2}")
                    lis_folders = os.popen("ls").readlines()
                    print(lis_folders)
                    n = 10 * len(lis_folders)
                    for i in range(len(lis_folders)):
                        os.chdir(f"{Saving_path}/{lis_folders[i][:-1]}")
                        lis_files = ["0_100", "100_200", "200_300", "300_400", "400_500", "500_600", "600_700",
                                     "700_800", "800_900", "900_1000"]

                        for j in range(10):
                            pm.cmd.load(f"./average{lis_files[j]}.pdb")
                            pm.cmd.color("red", "ss h")
                            pm.cmd.color("blue", "ss s")
                            pm.cmd.color("green", "ss l")
                            pm.cmd.bg_color(color='White')
                            pm.cmd.png(f"{pro_name1}_{lis_folders[i][:-1]}_{lis_files[j]}_{pro_name2}.png")
                            pm.cmd.delete("all")

                            os.system(f"mv {pro_name1}_{lis_folders[i][:-1]}_{lis_files[j]}_{pro_name2}.png {Saving_path}/.Pics_{pro_name1}_{pro_name2}")

                    ############################################################################################################
                    # Making a folder to moving the files and folders in that folder.
                    os.chdir(f"{source_path}/{source_folder_name}")     # Name
                    pro_name = os.popen("ls").readlines()
                    pro_name = pro_name[number_of_protein]
                    for i in range(len(pro_name)):
                        if pro_name[i] == "_":
                            pro_name1 = pro_name[0:i]
                            pro_name2 = pro_name[i + 1:-1]
                            break

                    os.chdir(Saving_path)
                    files = os.popen("ls -a").readlines()
                    files = files[2:]
                    for i in range(len(files)):
                        files[i] = files[i][:-1]

                    Analysis_folders = pro_name1 + "_10Frame_" + pro_name2
                    os.system(f"mkdir {Analysis_folders}")

                    for i in files:
                        print(f"Hello{i}")
                        print(f"{Analysis_folders}")
                        os.system(f"mv {i} {Analysis_folders}")

                    print(f"={Analysis_folders}")
                    print(type(Analysis_folders))
                    print(Analysis_folders == "1fsv_RMSF_0to1ns_2")
                    print(f"={final_Saving_path}")
                    print(type(final_Saving_path))
                    print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                    # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                    print(f"mv {Analysis_folders} {final_Saving_path}")
                    os.system(f"mv {Analysis_folders} {final_Saving_path}")

                    number_of_protein += 1

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style=3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style=7)} to exit the program) \n")

    elif alys_name == "hbond":                                                                                        # The Hydrogen bond numbers analysis
        ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'PLEASE READ', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        if usr == "y":
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            Saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if Saving_path == "def":
                Saving_path = f"/home/{nameOfUser}/Documents/Test3"

            print("We are now OPERATING on following path:", os.getcwd())  # The current path
            print(" ")

            source_path = input(f"What PATH do you in mind\
                                (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2"

            os.chdir(source_path)
            source_folder_name = os.popen("ls").read()[:-1]

            os.chdir(f"{source_path}/{source_folder_name}")
            proteins_path = os.popen("ls").readlines()
            for i in range(len(proteins_path)):
                proteins_path[i] = proteins_path[i][:-1]

            number_of_protein = 0
            for protein_in_path in proteins_path:
                ###########################################################################################################
                # The main work of Analysis and gromacs commands.
                # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                print("We are now OPERATING on following path:", os.getcwd())

                print(" ")  # Second folders name list
                print(" ")
                print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                os.popen(f'ls').readlines()

                # Finding Inputs(.tpr)
                files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_tpr)

                print(" ")
                print(" ")
                list_files_root_tpr.sort()
                print(list_files_root_tpr)
                print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                # Finding Inputs(.xtc)
                files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_xtc)

                print(" ")
                print(" ")
                list_files_root_xtc.sort()
                print(list_files_root_xtc)
                print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                print(" ")
                print(" ")
                c = 0
                index = []
                while c <= len(list_files_root_tpr):  # Making folders
                    os.chdir(Saving_path)
                    c += 1
                    if c == len(list_files_root_tpr) + 1:
                        break
                    os.system(f"mkdir hbond_{c}")
                    os.chdir(f"{Saving_path}/hbond_{c}")
                    index.append(f"hbond_{c}")
                print(index)

                print(" ")
                print(" ")
                each_root_tpr = files_roots_tpr.split()  # Copying the inputs to eac folders that we made(.tpr)
                print(each_root_tpr)
                length = len(each_root_tpr)
                i = 0
                list_tpr = []
                while i <= length - 1:
                    list_roots = each_root_tpr[i]  # .tpr
                    print(list_roots[1:])
                    print(" ")
                    print(" ")
                    each_folder = each_root_tpr[i].split("/")
                    list_tpr.append(each_folder[-1])
                    print(each_folder)
                    print(" ")
                    print(" ")
                    print(len(each_folder[-1]))
                    print(" ")
                    print(" ")
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    i += 1

                print(" ")
                print(" ")
                each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                print(each_root_xtc)
                length2 = len(each_root_xtc)
                i = 0
                list_xtc = []

                The_end_of_the_bar = (len(list_files_root_tpr)) * 5
                counting = 0

                while i <= length2 - 1:
                    list_roots = each_root_xtc[i]  # .xtc
                    print(list_roots[1:])
                    each_folder = each_root_xtc[i].split("/")
                    list_xtc.append(each_folder[-1])
                    print(each_folder)
                    print(len(each_folder[-1]))
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    os.chdir(Saving_path + "/" + index[i])

                    gmx = gmxExe
                    os.system(
                        f"echo 1 1 | {gmx} hbond -f {list_xtc[i]} -s {list_tpr[i]} -num hbond.xvg")  # Operating

                    counting += 1

                    The_Processor(len(list_files_root_tpr), counting)

                    hbond = open(f"hbond.xvg").readlines()
                    length12 = len(hbond)

                    time_list = []
                    hbond_list = []
                    index_1 = []
                    c = 0

                    for k in range(25, length12):
                        if hbond[k][8] == " ":
                            time_list.append(float(hbond[k][9]))
                        elif hbond[k][7] == " ":
                            time_list.append(float(hbond[k][8:10]))
                        elif hbond[k][6] == " ":
                            time_list.append(float(hbond[k][7:10]))
                        else:
                            time_list.append(float(hbond[k][6:10]))

                        if hbond[k][20] == " ":
                            hbond_list.append(float(hbond[k][21]))
                        elif hbond[k][19] == " ":
                            hbond_list.append(float(hbond[k][20:22]))
                        else:
                            hbond_list.append(float(hbond[k][19:22]))

                        c += 1
                        index_1.append(c)

                    ex_dict = {"time": time_list, "H_bond": hbond_list}

                    n = pd.DataFrame(data=ex_dict, index=index_1)

                    writer = pd.ExcelWriter('The_analysis_H_bond.xlsx', engine='xlsxwriter')
                    n.to_excel(writer, sheet_name="The H_bond Analysis")
                    writer.save()

                    i += 1
                    print(" ")
                    print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                ############################################################################################################
                # Changing the name of the folders of each Acceleration simulation's analysis
                ans = Saving_path
                dir = ans
                os.chdir(dir)
                folders = os.popen("ls").readlines()

                print(folders)

                for i in folders:
                    os.chdir(f"{dir}/hbond_{i[6:-1]}")
                    file = os.popen("ls").readlines()
                    print(file)
                    acc_num = file[1][7:-5]  # the first number is the number of files in each folders.
                    print(acc_num + "<")
                    if acc_num == "":
                        acc_num = 0
                    os.chdir(f"{dir}")
                    os.system(f"mv hbond_{i[6:-1]} hbond_{acc_num}")

                ############################################################################################################
                # Making a excel folder for the whole Analysis of each Acceleration simulation and partitioning it.
                ans2 = Saving_path
                os.chdir(ans2)

                lis_folders = os.popen("ls").readlines()
                print(lis_folders)
                H0_1000 = []
                H800_1000 = []
                dic_datas0_1000 = {}
                dic_datas800_1000 = {}
                index = []
                for i in range(len(lis_folders)):

                    os.chdir(f"{ans2}/{lis_folders[i][:-1]}")
                    lis_files = os.popen("ls").readlines()
                    for j in range(len(lis_files)):
                        lis_files[j] = lis_files[j][:-1]
                    print(lis_files)

                    datas = pd.read_excel("The_analysis_H_bond.xlsx")
                    c0_1000 = 0
                    sum0_1000H = 0
                    while c0_1000 != 101:
                        sum0_1000H += datas["H_bond"][c0_1000]
                        c0_1000 += 1
                    ave0_1000H = sum0_1000H / 100

                    H0_1000.append(ave0_1000H)

                    c800_1000 = 80
                    sum800_1000H = 0
                    while c800_1000 != 101:
                        sum800_1000H += datas["H_bond"][c800_1000]
                        c800_1000 += 1
                    ave800_1000rg = sum800_1000H / 100

                    H800_1000.append(ave800_1000rg)
                    index.append(lis_folders[i][6:-1])

                dic_datas0_1000["H_bond"] = H0_1000
                dic_datas800_1000["H_bond"] = H800_1000
                print(dic_datas0_1000)
                print(dic_datas800_1000)

                os.chdir(ans2)

                dic0_1000 = pd.DataFrame(data=dic_datas0_1000, index=index)
                print(dic0_1000)
                writer1 = pd.ExcelWriter('Total_Hbond.xlsx', engine='xlsxwriter')
                dic0_1000.to_excel(writer1, sheet_name="The H_bond Analysis")
                writer1.save()

                dic800_1000 = pd.DataFrame(data=dic_datas800_1000, index=index)
                writer2 = pd.ExcelWriter('800_1000_Hbond.xlsx', engine='xlsxwriter')
                dic800_1000.to_excel(writer2, sheet_name="The H_bond Analysis")
                writer2.save()
                ############################################################################################################
                # Making a folder to moving the files and folders in that folder.
                os.chdir(f"{source_path}/{source_folder_name}")
                pro_name = os.popen("ls").readlines()
                pro_name = pro_name[number_of_protein]
                for i in range(len(pro_name)):
                    if pro_name[i] == "_":
                        pro_name1 = pro_name[0:i]
                        pro_name2 = pro_name[i + 1:-1]
                        break

                os.chdir(Saving_path)
                files = os.popen("ls").readlines()
                for i in range(len(files)):
                    files[i] = files[i][:-1]

                Analysis_folders = pro_name1 + "_HBOND_" + pro_name2
                os.system(f"mkdir {Analysis_folders}")

                for i in files:
                    print(f"Hello{i}")
                    print(f"{Analysis_folders}")
                    os.system(f"mv {i} {Analysis_folders}")

                print(f"={Analysis_folders}")
                print(type(Analysis_folders))
                print(Analysis_folders == "1fsv_RMSD_0to1ns_2")
                print(f"={final_Saving_path}")
                print(type(final_Saving_path))
                print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                print(f"mv {Analysis_folders} {final_Saving_path}")
                os.system(f"mv {Analysis_folders} {final_Saving_path}")

                number_of_protein += 1
                ############################################################################################################

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style=3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style=7)} to exit the program) \n")

    elif alys_name == "rg":                                                                                           # The the radius of gyration analysis
        ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'PLEASE READ', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        if usr == "y":
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            Saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if Saving_path == "def":
                Saving_path = f"/home/{nameOfUser}/Documents/Test3"

            print("We are now OPERATING on following path:", os.getcwd())  # The current path
            print(" ")

            source_path = input(f"What PATH do you in mind\
                    (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2"

            os.chdir(source_path)
            source_folder_name = os.popen("ls").read()[:-1]

            os.chdir(f"{source_path}/{source_folder_name}")
            proteins_path = os.popen("ls").readlines()
            for i in range(len(proteins_path)):
                proteins_path[i] = proteins_path[i][:-1]

            number_of_protein = 0
            for protein_in_path in proteins_path:
                ###########################################################################################################
                # The main work of Analysis and gromacs commands.
                # First it make the folders and put the .xtc and .tpr files in them and do the analysis. the analysis is for ["0_200", "200_400", "400_600", "600_800", "800_1000", "950_1000"] ps of simulation.
                os.chdir(f"{source_path}/{source_folder_name}/{protein_in_path}")
                print("We are now OPERATING on following path:", os.getcwd())

                print(" ")  # Second folders name list
                print(" ")
                print(f"The files in {source_path}/{source_folder_name}/{protein_in_path} are:")
                os.popen(f'ls').readlines()

                # Finding Inputs(.tpr)
                files_roots_tpr = os.popen(f'find -type f -name "md*.tpr"').read()
                list_files_root_tpr = os.popen(f'find -type f -name "md*.tpr"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_tpr)

                print(" ")
                print(" ")
                list_files_root_tpr.sort()
                print(list_files_root_tpr)
                print(f"The Number of files(Inputs): {len(list_files_root_tpr)}")  # Making a list of roots

                # Finding Inputs(.xtc)
                files_roots_xtc = os.popen(f'find -type f -name "md*.xtc"').read()
                list_files_root_xtc = os.popen(f'find -type f -name "md*.xtc"').readlines()
                print(" ")
                print(" ")
                print("The ROOT of your wanted files: ")
                print(files_roots_xtc)

                print(" ")
                print(" ")
                list_files_root_xtc.sort()
                print(list_files_root_xtc)
                print(f"The Number of files(Inputs): {len(list_files_root_xtc)}")  # Making a list of roots

                print(" ")
                print(" ")
                c = 0
                index = []
                while c <= len(list_files_root_tpr):  # Making folders
                    os.chdir(Saving_path)
                    c += 1
                    if c == len(list_files_root_tpr) + 1:
                        break
                    os.system(f"mkdir rg_{c}")
                    os.chdir(f"{Saving_path}/rg_{c}")
                    index.append(f"rg_{c}")
                print(index)

                print(" ")
                print(" ")
                each_root_tpr = files_roots_tpr.split()  # Copying the inputs to eac folders that we made(.tpr)
                print(each_root_tpr)
                length = len(each_root_tpr)
                i = 0
                list_tpr = []
                while i <= length - 1:
                    list_roots = each_root_tpr[i]  # .tpr
                    print(list_roots[1:])
                    print(" ")
                    print(" ")
                    each_folder = each_root_tpr[i].split("/")
                    list_tpr.append(each_folder[-1])
                    print(each_folder)
                    print(" ")
                    print(" ")
                    print(len(each_folder[-1]))
                    print(" ")
                    print(" ")
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    i += 1

                print(" ")
                print(" ")
                each_root_xtc = files_roots_xtc.split()  # Copying the inputs to eac folders that we made(.xtc)
                print(each_root_xtc)
                length2 = len(each_root_xtc)
                i = 0
                list_xtc = []

                The_end_of_the_bar = (len(list_files_root_tpr)) * 5
                counting = 0

                while i <= length2 - 1:
                    list_roots = each_root_xtc[i]  # .xtc
                    print(list_roots[1:])
                    each_folder = each_root_xtc[i].split("/")
                    list_xtc.append(each_folder[-1])
                    print(each_folder)
                    print(len(each_folder[-1]))
                    os.chdir(
                        f"{source_path}/{source_folder_name}/{protein_in_path}/{list_roots[1:-(len(each_folder[-1]))]}")
                    os.system(f"cp {each_folder[-1]} {Saving_path}/{index[i]}")
                    print(Saving_path + "/" + index[i])

                    os.chdir(Saving_path + "/" + index[i])

                    gmx = gmxExe
                    os.system(
                        f"echo 1 | {gmx} gyrate -f {list_xtc[i]} -s {list_tpr[i]} -o {list_xtc[i][7:-4]}.xvg")  # Operating

                    counting += 1

                    The_Processor(len(list_files_root_tpr), counting)

                    rg = open(f"{list_xtc[i][7:-4]}.xvg").readlines()
                    print(rg)
                    length12 = len(rg)

                    time_list = []
                    rg_list = []
                    rgx_list = []
                    rgy_list = []
                    rgz_list = []
                    index_1 = []
                    c = 0

                    for k in range(27, length12):
                        # time
                        print(rg[k][:9])
                        time_list.append(float(rg[k][:10]))

                        # rg
                        print(rg[k][15:21])
                        rg_list.append(float(rg[k][15:22]))

                        # rgx
                        print(rg[k][27:33])
                        rgx_list.append(float(rg[k][27:34]))

                        # rgy
                        print(rg[k][39:45])
                        rgy_list.append(float(rg[k][39:46]))

                        # rgz
                        print(rg[k][51:])
                        rgz_list.append(float(rg[k][51:]))

                        c += 1
                        index_1.append(c)

                    ex_dict = {"time": time_list, "Rg": rg_list, "Rg_x": rgx_list, "Rg_y": rgy_list, "Rg_z": rgz_list}

                    n = pd.DataFrame(data=ex_dict, index=index_1)

                    writer = pd.ExcelWriter('The_analysis_Rg.xlsx', engine='xlsxwriter')
                    n.to_excel(writer, sheet_name="The Rg Analysis")
                    writer.save()

                    print(i)
                    i += 1

                    print(" ")
                    print(f"{ansi(42, 'The Process is COMPLETE!!!', style=6)}")

                ############################################################################################################
                # Changing the name of the folders of each Acceleration simulation's analysis
                ans = Saving_path
                dir = ans
                os.chdir(dir)
                folders = os.popen("ls").readlines()

                print(folders)

                for i in folders:
                    os.chdir(f"{dir}/rg_{i[3:-1]}")
                    file = os.popen("ls").readlines()
                    print(file)
                    acc_num = file[2][7:-5]  # the first number is the number of files in each folders.
                    print(acc_num + "<")
                    if acc_num == "":
                        acc_num = 0
                    os.chdir(f"{dir}")
                    os.system(f"mv rg_{i[3:-1]} rg_{acc_num}")

                ############################################################################################################
                # Making a excel folder for the whole Analysis of each Acceleration simulation and partitioning it.
                ans2 = Saving_path
                os.chdir(ans2)

                lis_folders = os.popen("ls").readlines()
                print(lis_folders)
                Rg0_1000 = []
                Rgx0_1000 = []
                Rgy0_1000 = []
                Rgz0_1000 = []
                Rg800_1000 = []
                Rgx800_1000 = []
                Rgy800_1000 = []
                Rgz800_1000 = []
                dic_datas0_1000 = {}
                dic_datas800_1000 = {}
                index = []
                for i in range(len(lis_folders)):

                    os.chdir(f"{ans2}/{lis_folders[i][:-1]}")
                    lis_files = os.popen("ls").readlines()
                    for j in range(len(lis_files)):
                        lis_files[j] = lis_files[j][:-1]
                    print(lis_files)

                    datas = pd.read_excel("The_analysis_Rg.xlsx")
                    c0_1000 = 0
                    sum0_1000rg = 0
                    sum0_1000rgx = 0
                    sum0_1000rgy = 0
                    sum0_1000rgz = 0
                    while c0_1000 != 101:
                        sum0_1000rg += datas["Rg"][c0_1000]
                        sum0_1000rgx += datas["Rg_x"][c0_1000]
                        sum0_1000rgy += datas["Rg_y"][c0_1000]
                        sum0_1000rgz += datas["Rg_z"][c0_1000]
                        c0_1000 += 1
                    ave0_1000rg = sum0_1000rg / 100
                    ave0_1000rgx = sum0_1000rgx / 100
                    ave0_1000rgy = sum0_1000rgy / 100
                    ave0_1000rgz = sum0_1000rgz / 100

                    Rg0_1000.append(ave0_1000rg)
                    Rgx0_1000.append(ave0_1000rgx)
                    Rgy0_1000.append(ave0_1000rgy)
                    Rgz0_1000.append(ave0_1000rgz)

                    c800_1000 = 80
                    sum800_1000rg = 0
                    sum800_1000rgx = 0
                    sum800_1000rgy = 0
                    sum800_1000rgz = 0
                    while c800_1000 != 101:
                        sum800_1000rg += datas["Rg"][c800_1000]
                        sum800_1000rgx += datas["Rg_x"][c800_1000]
                        sum800_1000rgy += datas["Rg_y"][c800_1000]
                        sum800_1000rgz += datas["Rg_z"][c800_1000]
                        c800_1000 += 1
                    ave800_1000rg = sum800_1000rg / 100
                    ave800_1000rgx = sum800_1000rgx / 100
                    ave800_1000rgy = sum800_1000rgy / 100
                    ave800_1000rgz = sum800_1000rgz / 100

                    Rg800_1000.append(ave800_1000rg)
                    Rgx800_1000.append(ave800_1000rgx)
                    Rgy800_1000.append(ave800_1000rgy)
                    Rgz800_1000.append(ave800_1000rgz)

                    index.append(lis_folders[i][3:-1])

                dic_datas0_1000["Rg"] = Rg0_1000
                dic_datas0_1000["Rgx"] = Rgx0_1000
                dic_datas0_1000["Rgy"] = Rgy0_1000
                dic_datas0_1000["Rgz"] = Rgz0_1000
                dic_datas800_1000["Rg"] = Rg800_1000
                dic_datas800_1000["Rgx"] = Rgx800_1000
                dic_datas800_1000["Rgy"] = Rgy800_1000
                dic_datas800_1000["Rgz"] = Rgz800_1000

                print(dic_datas0_1000)
                print(dic_datas800_1000)

                os.chdir(ans2)
                dic0_1000 = pd.DataFrame(data=dic_datas0_1000, index=index)
                print(dic0_1000)
                writer1 = pd.ExcelWriter('Total_Rg.xlsx', engine='xlsxwriter')
                dic0_1000.to_excel(writer1, sheet_name="The Rg Analysis")
                writer1.save()

                dic800_1000 = pd.DataFrame(data=dic_datas800_1000, index=index)
                writer2 = pd.ExcelWriter('800_1000_Rg.xlsx', engine='xlsxwriter')
                dic800_1000.to_excel(writer2, sheet_name="The Rg Analysis")
                writer2.save()
                ############################################################################################################
                # Making a folder to moving the files and folders in that folder.
                os.chdir(f"{source_path}/{source_folder_name}")
                pro_name = os.popen("ls").readlines()
                pro_name = pro_name[number_of_protein]
                for i in range(len(pro_name)):
                    if pro_name[i] == "_":
                        pro_name1 = pro_name[0:i]
                        pro_name2 = pro_name[i + 1:-1]
                        break

                os.chdir(Saving_path)
                files = os.popen("ls").readlines()
                for i in range(len(files)):
                    files[i] = files[i][:-1]

                Analysis_folders = pro_name1 + "_Rg_" + pro_name2
                os.system(f"mkdir {Analysis_folders}")

                for i in files:
                    print(f"Hello{i}")
                    print(f"{Analysis_folders}")
                    os.system(f"mv {i} {Analysis_folders}")

                print(f"={Analysis_folders}")
                print(type(Analysis_folders))
                print(Analysis_folders == "1fsv_RMSD_0to1ns_2")
                print(f"={final_Saving_path}")
                print(type(final_Saving_path))
                print(final_Saving_path == f"/home/{nameOfUser}/Documents/Test1")
                # os.system(f"sudo chmod +x {Saving_path}/{Analysis_folders}")
                print(f"mv {Analysis_folders} {final_Saving_path}")
                os.system(f"mv {Analysis_folders} {final_Saving_path}")

                number_of_protein += 1
                ############################################################################################################

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style=3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style=7)} to exit the program) \n")
    # elif alys_name == "trjconv":                                                                                      # The Trajectory analysis

    elif alys_name == "cm":
        ############################################################################################################
        # Inputing the paths for saving and working
        print(f"                  {ansi(34, 'PLEASE READ', style=4)}")
        print("Make sure doing these work before using the app:")
        print(f"1. Set the {ansi(4, '.xtc')} and {ansi(4, '.tpr')} and {ansi(4, '.gro')} files being seperated")
        print("   from DOING_Directory, exactly three folder away.")
        print("2. You should doing it in your computer directories and not in hard driver or sth.")

        usr = input("Do you want to continue? (y / 'Enter')")

        if usr == "y":
            # Working and Saving Points and Directories
            final_Saving_path = input(
                f"Where do you want to save your all your outputs?(ex: '/home/{nameOfUser}/Documents/Test1' = def): ")
            if final_Saving_path == "def":
                final_Saving_path = f"/home/{nameOfUser}/Documents/Test1"

            saving_path = input(
                f"Where do you want to SAVE your data? (Output)(ex: '/home/{nameOfUser}/Documents/Test3' = def): ")  # The saving path
            if saving_path == "def":
                saving_path = f"/home/{nameOfUser}/Documents/Test3"

            source_path = input(f"What PATH do you in mind\
                    (YOUR INPUTS)(Proteins-folder's path)(ex: '/home/{nameOfUser}/Documents/Test2/proteins' = def)?: ")  # The source of Inputs
            if source_path == "def":
                source_path = f"/home/{nameOfUser}/Documents/Test2/proteins"

            # Making the list of Acceleration for analysis
            listAccelerations = input("Which Accelerations do you want to analysis?(first pick the type of formatting)('list' or 'num')\n"
                                      "\tThe Instruction to write it:\n"
                                      "\t\tIf you pick : 'list' ---example---> > 0.00-0.30 0.28 0.03-0.21\n"
                                      "\t\tIf you pick : 'num' ---example----> > 0.00-0.30 or > 0.28\n")
            if listAccelerations == "list":
                listAcc = input("Now enter the list of Accelerations:\n> ")
                listAccelerations = listAcc.split(" ")

            elif listAccelerations == "num":
                listAccelerations = []
                numAcc = input("Now enter the Accelerations one by one:\n> ")
                while numAcc != "":
                    listAccelerations.append(numAcc)
                    numAcc = input("Enter the next Acceleration:(if you want add no more, press Enter)\n> ")

            # Setting the beggining frame and final frame:
            bFrame = int(input("Enter the Beginning Frame (Usually 0)\n> "))
            fFrame = int(input("Enter the Final Frame (Usually 100)\n> "))

            # Setting the Cmap Style:
            cmapStyle = input("What style do you want?(gnuplot, seismic, 'Accent', 'Accent_r', 'Blues', 'Blues_r',\n"
                              " 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2',\n"
                              " 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r',\n"
                              " 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r',\n"
                              " 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr',\n"
                              " 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r',\n"
                              " 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1',\n"
                              " 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r',\n"
                              " 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',\n"
                              " 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r',\n"
                              " 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper',\n"
                              " 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r',\n"
                              " 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',\n"
                              " 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r',\n"
                              " 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv',\n"
                              " 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral',\n"
                              " 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism',\n"
                              " 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer',\n"
                              " 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r',\n"
                              " 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted',\n"
                              " 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r')\n> ")

            # Setting the Cutoff of the analysis:
            cutoff = float(input("What value (in Angstrom) for CUTOFF do you want to set? (For Hydrogen bonds : 1.5 Angestrom)\n> "))

            # The Analysis:
            os.chdir(f"{source_path}")
            simulationFolders = listOfFilesAndFolders()
            print(simulationFolders)

            c = 0
            n = len(simulationFolders) * 2 + len(listAccelerations) * 3 * len(simulationFolders)
            The_Processor(n, c)

            for simulationFolder in simulationFolders:
                c += 1
                The_Processor(n, c)
                print(f"!!!!!!!!!!!!{simulationFolder}!!!!!!!!!!!!")
                os.chdir(f"{source_path}/{simulationFolder}")
                print(listOfFilesAndFolders())
                simulationAccFolders = listOfFilesAndFolders()
                for acc in listAccelerations:
                    os.chdir(f"{source_path}/{simulationFolder}")
                    print(f"!!!!!!!!!!!!{acc}!!!!!!!!!!!!")
                    os.system(f"mkdir {saving_path}/{acc}")
                    accelerations = acc.split("-")
                    print(f"!!!!!!!!!!!!{accelerations}!!!!!!!!!!!!")
                    typeOfContactMap = len(accelerations)
                    xtcFilesDir = []
                    tprFilesDir = []
                    groFilesDir = []
                    c += 1
                    The_Processor(n, c)
                    for file in accelerations:
                        print(file)
                        print(os.popen(f'find -type f -name "md*{file}.xtc"').read()[:-1])
                        xtcFilesDir.append(os.popen(f'find -type f -name "md*{file}.xtc"').read()[:-1])
                        tprFilesDir.append(os.popen(f'find -type f -name "md*{file}.tpr"').read()[:-1])
                        groFilesDir.append(os.popen(f'find -type f -name "md*{file}.gro"').read()[:-1])
                    print(xtcFilesDir)
                    print(tprFilesDir)
                    print(groFilesDir)
                    for ind in range(len(xtcFilesDir)):
                        print(
                            f"{xtcFilesDir[ind][:2]}\\{xtcFilesDir[ind][2:7]}\\{xtcFilesDir[ind][7:]} {saving_path}/{acc}")
                        copyFiles(f"{xtcFilesDir[ind][:2]}\\{xtcFilesDir[ind][2:7]}\\{xtcFilesDir[ind][7:]}",
                                  f"{saving_path}/{acc}")
                        copyFiles(f"{tprFilesDir[ind][:2]}\\{tprFilesDir[ind][2:7]}\\{tprFilesDir[ind][7:]}",
                                  f"{saving_path}/{acc}")
                        copyFiles(f"{groFilesDir[ind][:2]}\\{groFilesDir[ind][2:7]}\\{groFilesDir[ind][7:]}",
                                  f"{saving_path}/{acc}")

                    os.chdir(f"{saving_path}/{acc}")
                    if typeOfContactMap == 1:
                        c += 1
                        The_Processor(n, c)
                        os.system(
                            f"echo 1 | {gmxExe} trjconv -f md.acc_{accelerations[0]}.xtc -s md.acc_{accelerations[0]}.tpr -o protein_{accelerations[0]}.xtc")
                        os.system(
                            f"echo 1 | {gmxExe} trjconv -f md.acc_{accelerations[0]}.gro -s md.acc_{accelerations[0]}.tpr -o protein_{accelerations[0]}.gro")
                        c += 1
                        The_Processor(n, c)

                        traj = md.load(f"protein_{accelerations[0]}.xtc", top=f"protein_{accelerations[0]}.gro")
                        topology = traj.topology
                        # switch1 = topology.select("resid 0 to 577")

                        # frame_contact = ContactFrequency(traj[0], query=switch1)
                        print(traj[bFrame:fFrame])
                        frame_contact = ContactFrequency(traj[bFrame:fFrame],
                                                         cutoff=cutoff)  # traj[n] is about the nth frame

                        # frame_contact.save_to_file(f"traj_contacts_{num}.p")

                        (fig, ax) = frame_contact.residue_contacts.plot(cmap=cmapStyle, vmin=-1.5, vmax=1.5)
                        plt.xlabel("Residue")
                        plt.ylabel("Residue")

                        plt.savefig(f"{simulationFolder[:4]}_{simulationFolder[4:]}_contacts_{accelerations[0]}.png")

                        d = frame_contact.residue_contacts.df
                        dd = pd.DataFrame(d)
                        dd.to_excel(f"{simulationFolder[:4]}_{simulationFolder[4:]}_contacts_{accelerations[0]}.xlsx")

                    elif typeOfContactMap == 2:     # The Difference between 'first Acc' to  'second Acc' ---> 'firstAcc-secondAcc'
                        traj = []
                        map = []
                        for i in accelerations:
                            c += 1
                            The_Processor(n, c)
                            os.system(f"echo 1 | {gmxExe} trjconv -f md.acc_{i}.xtc -s md.acc_{i}.tpr -o protein_{i}.xtc")
                            os.system(f"echo 1 | {gmxExe} trjconv -f md.acc_{i}.gro -s md.acc_{i}.tpr -o protein_{i}.gro")
                            The_Processor(n, c)
                            t = md.load(f"protein_{i}.xtc", top=f"protein_{i}.gro")
                            traj.append(t)
                            m = ContactFrequency(t[bFrame:fFrame], cutoff=cutoff)
                            map.append(m)
                        print(map)
                        diff = otcd(map[0], map[1], topology=traj[0].topology)                      # like this: >>> diff = map[0] - map[1] 
                        (fig, ax) = diff.residue_contacts.plot(cmap=cmapStyle, vmin=-1.5, vmax=1.5) #            >>> (fig, ax) = diff.residue_... 

                        plt.xlabel("Residue")
                        plt.ylabel("Residue")

                        fig.savefig(f"{simulationFolder[:4]}_{simulationFolder[4:]}_DiffContacts_{acc}.png")

                # The Folder making for each simulation
                os.chdir(final_Saving_path)
                folderName = f"{simulationFolder[:4]}_ContactMap_{simulationFolder[4:]}"
                os.system(f"mkdir {folderName}")

                # The Moving the files from 'saving directory' to 'final saving directory'
                os.chdir(saving_path)
                files = listOfFilesAndFolders()
                for file in files:
                    moveFiles(file, f"{final_Saving_path}/{folderName}")
                c += 1
                The_Processor(n, c)

        elif usr != "n":
            alys_name = input(
                f"Do you want do another {ansi(1, 'analysis')}?\n(The {ansi(5, 'options', style=3)}: '{ansi(31, 'dssp')}', '{ansi(32, 'rmsd')}', '{ansi(33, 'rmsf')}', '{ansi(34, 'hbond')}', '{ansi(35, 'rg')}', '{ansi(36, 'trjconv')}', '{ansi(30, 'cm')}'(Contact Map), or type just {ansi(94, 'Enter', style=7)} to exit the program) \n")

    elif alys_name == "":
        print(f"{ansi(33, 'See you next time!!!', style = 1)}")

    else:
        print(f"{ansi(31, 'Your name choice is not in the options!', style = 4)}")
