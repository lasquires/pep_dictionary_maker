import os.path
import pprint
#this is a list of all amino acids and their neutral masses.
global aa_mass
aa_mass = {'A': 71.037114, 'R':156.101111 ,'N':	114.042927, 'D':	115.026943, 'C':	103.009185, 'E':	129.042593, 'Q'	: 128.058578,\
              'G':	57.021464, 'H':	137.058912, 'I' :	113.084064, 'L':	113.084064, 'K':	128.094963, 'M' :	131.040485, 'F':	147.068414,\
              'P':	97.052764, 'S':	87.032028, 'T':	101.047679, 'U':	150.95363, 'W':	186.079313, 'Y':	163.06332, 'V':	99.068414}
global peptideDict
peptideDict = {} #I'm sorry, I know you all hate global var, but it seemed necessary here...


def File_to_Dict(filename): #this lets you read in an already existing dictionary:
    if not os.path.isfile(filename):
        return {}
    with open(filename) as f:
        d = {}
        for line in f:
            if line[0].isdigit():
                key, value = line.strip().split(': ')
                d[key] = [i.strip("'") for i in value.strip('[]').split(', ')]
        return d

def to_File(dict, filename): #this lets you write the dictionary to the file
    filename += '.txt'
    with open(filename, 'w') as f:
        for key in dict:
            f.write(str(key) + ': ' + str(dict[key]) + '\n')

def reverse_dict():
    reversed = {}
    for k in peptideDict:
        if peptideDict[k] not in reversed.keys():
            reversed[peptideDict[k]] =[]
        reversed[peptideDict[k]].append(k)
    myKeys = list(reversed.keys())
    myKeys.sort()
    sorted_dict = {i: reversed[i] for i in myKeys}
    reversed = sorted_dict
    return reversed

def sort_dict(dict):
    myKeys = list(dict.keys())
    floatedKeys = []
    for i in myKeys:
        floatedKeys.append(float(i))
    floatedKeys.sort()
    #myKeys.sort()
    sorted_dict = {str(i): dict[str(i)] for i in floatedKeys}
    return sorted_dict

def generate_combinations(length, prefix):
    if len(prefix) == length:
        mass_Hydrogen = 1.0078
        mass_Oxygen = 15.994915
        mass_proton = 1.007276466879
        m_z_1 = mass_proton + mass_Oxygen + 2 * mass_Hydrogen
        for aa in prefix:
            m_z_1 += aa_mass[aa]
        peptideDict[prefix] = round(m_z_1, 4)
    else:
        for char in aa_mass:
            #if char not in prefix:
                # check if the prefix is not empty and the last character of prefix is greater than the current char
                # this condition helps to eliminate duplicates like 'abdc' and 'dbca'
                if prefix and prefix[-1] > char:
                    continue
                generate_combinations(length, prefix + char)

def add_to_dict(master_dict, length):
    peptideDict = {}
    generate_combinations(length, '')
    reversed = reverse_dict()
    for key in reversed:
        if str(key) not in master_dict.keys():
            master_dict[str(key)] = reversed[key]
        else:
            for peptide in reversed[key]:
                if peptide not in master_dict[str(key)]:
                    master_dict[str(key)].append(peptide)
    return sort_dict(master_dict)

def filter_by_presence(peptide, master_dict):
    filter_dict = {}
    for key in ultimate_dict:
        for pep in ultimate_dict[key]:
            canAdd = True
            for aa in pep:
                if aa not in peptide:
                    canAdd = False
        if canAdd:
            if key not in filter_dict.keys():
                filter_dict[key] = []
            filter_dict[key].append(pep)
    return filter_dict

def filter_by_m_z(mz_val, acceptance, master_dict):
    filter_dict = {}
    for key in ultimate_dict:
        key_as_num = float(key)
        if(abs(key_as_num - mz_val) <= acceptance):
            filter_dict[key] = ultimate_dict[key]
    return filter_dict

ultimate_dict = File_to_Dict('peptide_mz_dict.txt')
#
# peptide = ('MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFS\
# QYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMA\
# DCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYF\
# YAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGER\
# ALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQD\
# TISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGS\
# FLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIK\
# QNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCT\
# EDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTF\
# HADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFA\
# VEGPKLVVSTQTALA')
#
# filteredDict = filter_by_presence(peptide, dict)
# filteredDict = filter_by_m_z(353, 1, filteredDict)
# pprint.pprint(filteredDict)
#to_File(filteredDict, 'filtered')
max_peptide_combos = 9
for i in range(1,max_peptide_combos,1):
    ultimate_dict = add_to_dict(ultimate_dict, i+1)
print('This file has 100% of all peptide combinations of up to an mz value of ' + str(round((19.017791466879 + max_peptide_combos* aa_mass['G']), 4)))
to_File(ultimate_dict, 'peptide_mz_dict')





