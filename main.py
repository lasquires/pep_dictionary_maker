import os.path

aa_mass = {'A': 71.037114, 'R':156.101111 ,'N':	114.042927, 'D':	115.026943, 'C':	103.009185, 'E':	129.042593, 'Q'	: 128.058578,\
              'G':	57.021464, 'H':	137.058912, 'I' :	113.084064, 'L':	113.084064, 'K':	128.094963, 'M' :	131.040485, 'F':	147.068414,\
              'P':	97.052764, 'S':	87.032028, 'T':	101.047679, 'U':	150.95363, 'W':	186.079313, 'Y':	163.06332, 'V':	99.068414}
def check_mz(peptide):
    neutral_mass = 0
    for aa in peptide:
        neutral_mass += aa_mass[aa]
    mass_Hydrogen = 1.0078
    mass_Oxygen = 15.994915
    mass_proton = 1.007276466879
    mz = neutral_mass + 2*mass_Hydrogen +mass_Oxygen+mass_proton
    return mz
def File_to_Dict(filename):
    if not os.path.isfile(filename):
        return {}
    with open(filename) as f:
        d = {}
        for line in f:
            if line[0].isdigit():
                key, value = line.strip().split(': ')
                d[key] = [i.strip("'") for i in value.strip('[]').split(', ')]
        return d

def to_File(dict, filename):
  filename += '.txt'
  with open (filename , 'w') as f:
    for key in dict:
        f.write(str(key) + ': '+ str(dict[key]) + '\n')

def reverse_dict(dict):
    reversed = {}
    for k in dict:
        if dict[k] not in reversed.keys():
            reversed[dict[k]] =[]
        reversed[dict[k]].append(k)
    myKeys = list(reversed.keys())
    myKeys.sort()
    sorted_dict = {i: reversed[i] for i in myKeys}
    reversed = sorted_dict
    return reversed

def get_dimer(dict):
    dimer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            sorted_aa = ''.join(sorted([first,second]))
            mass_Hydrogen = 1.0078
            mass_Oxygen = 15.994915
            mass_proton = 1.007276466879
            neutral_mass = aa_mass[first] + aa_mass[second] + (2*mass_Hydrogen) + mass_Oxygen
            dimer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    dimer_dict = reverse_dict(dimer_dict)
    return dimer_dict

def get_trimer(dict):
    trimer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            for third in aa_mass:
                sorted_aa = ''.join(sorted([first,second,third]))
                mass_Hydrogen = 1.0078
                mass_Oxygen = 15.994915
                mass_proton = 1.007276466879
                neutral_mass = aa_mass[first] + aa_mass[second] + aa_mass[third] + (2*mass_Hydrogen) + mass_Oxygen
                trimer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    trimer_dict = reverse_dict(trimer_dict)
    return trimer_dict


def get_tetramer(dict):
    tetramer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            for third in aa_mass:
                for fourth in aa_mass:
                    sorted_aa = ''.join(sorted([first,second,third, fourth]))
                    mass_Hydrogen = 1.0078
                    mass_Oxygen = 15.994915
                    mass_proton = 1.007276466879
                    neutral_mass = aa_mass[first] + aa_mass[second] + aa_mass[third] + aa_mass[fourth] + (2 * mass_Hydrogen) + mass_Oxygen
                    tetramer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    tetramer_dict = reverse_dict(tetramer_dict)
    return tetramer_dict


def get_pentamer(dict):
    pentamer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            for third in aa_mass:
                for fourth in aa_mass:
                    for fifth in aa_mass:
                        sorted_aa = ''.join(sorted([first,second,third, fourth, fifth]))
                        mass_Hydrogen = 1.0078
                        mass_Oxygen = 15.994915
                        mass_proton = 1.007276466879
                        neutral_mass = aa_mass[first] + aa_mass[second] + aa_mass[third] + aa_mass[fourth] + aa_mass[fifth] + (2 * mass_Hydrogen) + mass_Oxygen
                        pentamer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    pentamer_dict = reverse_dict(pentamer_dict)
    return pentamer_dict


def get_hexamer(dict):
    hexamer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            for third in aa_mass:
                for fourth in aa_mass:
                    for fifth in aa_mass:
                        for sixth in aa_mass:
                            sorted_aa = ''.join(sorted([first,second,third, fourth, fifth, sixth]))
                            mass_Hydrogen = 1.0078
                            mass_Oxygen = 15.994915
                            mass_proton = 1.007276466879
                            neutral_mass = aa_mass[first] + aa_mass[second] + aa_mass[third] + aa_mass[fourth] + aa_mass[fifth] + aa_mass[sixth]+ (2 * mass_Hydrogen) + mass_Oxygen
                            hexamer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    hexamer_dict = reverse_dict(hexamer_dict)
    return hexamer_dict

def get_heptamer(dict):
    heptamer_dict = {}
    for first in aa_mass:
        for second in aa_mass:
            for third in aa_mass:
                for fourth in aa_mass:
                    for fifth in aa_mass:
                        for sixth in aa_mass:
                            for seventh in aa_mass:
                                sorted_aa = ''.join(sorted([first,second,third, fourth, fifth, sixth, seventh]))
                                mass_Hydrogen = 1.0078
                                mass_Oxygen = 15.994915
                                mass_proton = 1.007276466879
                                neutral_mass = aa_mass[first] + aa_mass[second] + aa_mass[third] + aa_mass[fourth] + aa_mass[fifth] + aa_mass[sixth]+ aa_mass[seventh]+(2 * mass_Hydrogen) + mass_Oxygen
                                heptamer_dict[sorted_aa] = round((neutral_mass + mass_proton), 4)
    heptamer_dict = reverse_dict(heptamer_dict)
    return heptamer_dict



def add_to_dict(mainDict, dictToAdd):
    for key in dictToAdd:
        if key not in mainDict.keys():
            mainDict[key] = dictToAdd[key]
        else:
            for peptide in dictToAdd[key]:
                mainDict[key].append(peptide)

##you don't have to run all at once.
##And yeah, sorry, recursion hurts my brain. copy and pasting is easier for just a few functions
ultimate_dict = File_to_Dict('peptide_dict.txt')
# add_to_dict(ultimate_dict, get_dimer(aa_mass))
# add_to_dict(ultimate_dict, get_trimer(aa_mass))
# add_to_dict(ultimate_dict, get_tetramer(aa_mass))
# add_to_dict(ultimate_dict, get_pentamer(aa_mass))
#add_to_dict(ultimate_dict, get_heptamer(aa_mass)) #MAKES THE FILE 20GB OR SOMETHING...PLZ DON'T DO
to_File(ultimate_dict, 'peptide_dict')





