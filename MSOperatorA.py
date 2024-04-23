from MSLogging import CLogging

from MSDataPeptide import CInfo1, CPeptideInfo, CAddModInfo, CModSite, CPeptideIndexSortInfo


def op_FillCInfo1(input_CEnzymeInfo:CInfo1, enzyme_info_str:str, length_up:int, length_low:int, mass_up:float, mass_low:float, digest_type:int=0, miss_clv:int=0):
    enzyme_info_list = enzyme_info_str.split(" ")
    flag = enzyme_info_list[2]
    if flag != "N" and flag != "C":
        print("[Error] Enzyme split flag can only be N or C.")
        info = "[Error] Enzyme split flag can only be N or C."
        CLogging.LogError(info)
    for i, aa in enumerate(enzyme_info_list[1]):
        input_CEnzymeInfo.enzyme_aa.append(aa)
        input_CEnzymeInfo.enzyme_flag.append(flag)

    if digest_type < 0 or digest_type> 2:
        print("[Error] Enzyme type can only be 0, 1, 2.")
        info = "[Error] Enzyme type can only be 0, 1, 2."
        CLogging.LogError(info)
    input_CEnzymeInfo.digest_type = digest_type
    if miss_clv < 0:
        print("[Error] Enzyme miss cleavage number must be >= 0.")
        info = "[Error] Enzyme miss cleavage number must be >= 0."
        CLogging.LogError(info)
    input_CEnzymeInfo.miss_clv = miss_clv

    if length_up < 0:
        info = "Peptide length maximum must be >= 0."
        CLogging.LogError(info)
    input_CEnzymeInfo.length_up = length_up

    if length_low < 0:
        info = "Peptide length minimum must be >= 0."
        CLogging.LogError(info)
    input_CEnzymeInfo.length_low = length_low

    if mass_up < 0:
        info = "Peptide mass maximum must be >= 0."
        CLogging.LogError(info)
    input_CEnzymeInfo.mass_up = mass_up

    if mass_low < 0:
        info = "Peptide mass minimum must be >= 0."
        CLogging.LogError(info)
    input_CEnzymeInfo.mass_low = mass_low

def op_UpdateCPeptideMathInfo(input_CPeptideInfo:CPeptideInfo, mass:float=0.0, gdm:float=0.0):
    input_CPeptideInfo.mass = mass
    input_CPeptideInfo.gdm = gdm

def op_FillCAddModInfo(input_CAddModInfo:CAddModInfo, fix_mod_name_list, ooput593_name_list, INI_DIC_MOD):
    m2mod = INI_DIC_MOD
    for mod_name in fix_mod_name_list:
        if len(mod_name) == 0:
            continue
        if mod_name not in m2mod:
            print("[Error] {} not in modification.ini".format(mod_name))
            return -1
        if mod_name in input_CAddModInfo.modname2index: return -1
        input_CAddModInfo.modname2index[mod_name] = len(input_CAddModInfo.mod_list)
        input_CAddModInfo.mod_list.append(m2mod[mod_name])
        input_CAddModInfo.modname2type[mod_name] = 0
        mod = m2mod[mod_name]
        if mod.type == "NORMAL":
            for a in mod.aa:
                if a in input_CAddModInfo.fix_aa: return -1
                input_CAddModInfo.fix_aa[a] = mod
        elif mod.type == "PEP_N":
            for a in mod.aa:
                if a in input_CAddModInfo.fix_N_aa: return -1
                input_CAddModInfo.fix_N_aa[a] = mod
        elif mod.type == "PEP_C":
            for a in mod.aa:
                if a in input_CAddModInfo.fix_C_aa: return -1
                input_CAddModInfo.fix_C_aa[a] = mod
        elif mod.type == "PRO_N":
            for a in mod.aa:
                if a in input_CAddModInfo.fix_PN_aa: return -1
                input_CAddModInfo.fix_PN_aa[a] = mod
        elif mod.type == "PRO_C":
            for a in mod.aa:
                if a in input_CAddModInfo.fix_PC_aa: return -1
                input_CAddModInfo.fix_PC_aa[a] = mod

    for mod_name in ooput593_name_list:
        if len(mod_name) == 0:
            continue
        if mod_name not in m2mod:
            print("[Error] {} not in modification.ini".format(mod_name))
            return -1
        if mod_name in input_CAddModInfo.modname2index: return -1
        input_CAddModInfo.modname2index[mod_name] = len(input_CAddModInfo.mod_list)
        input_CAddModInfo.mod_list.append(m2mod[mod_name])
        input_CAddModInfo.modname2type[mod_name] = 1
        mod = m2mod[mod_name]
        if mod.type == "NORMAL":
            for a in mod.aa:
                if a not in input_CAddModInfo.var_aa: input_CAddModInfo.var_aa[a] = []
                input_CAddModInfo.var_aa[a].append(mod)
        elif mod.type == "PEP_N":
            for a in mod.aa:
                if a not in input_CAddModInfo.var_N_aa: input_CAddModInfo.var_N_aa[a] = []
                input_CAddModInfo.var_N_aa[a].append(mod)
        elif mod.type == "PEP_C":
            for a in mod.aa:
                if a not in input_CAddModInfo.var_C_aa: input_CAddModInfo.var_C_aa[a] = []
                input_CAddModInfo.var_C_aa[a].append(mod)
        elif mod.type == "PRO_N":
            for a in mod.aa:
                if a not in input_CAddModInfo.var_PN_aa: input_CAddModInfo.var_PN_aa[a] = []
                input_CAddModInfo.var_PN_aa[a].append(mod)
        elif mod.type == "PRO_C":
            for a in mod.aa:
                if a not in input_CAddModInfo.var_PC_aa: input_CAddModInfo.var_PC_aa[a] = []
                input_CAddModInfo.var_PC_aa[a].append(mod)

def op_AddModIntoCPeptideInfo(input_CPeptideInfo:CPeptideInfo, mod:CModSite):

    input_CPeptideInfo.mods.append(mod)

def op_FillCPeptideIndexSortInfo(input_Data:CPeptideIndexSortInfo, pep_idx_folder:str, pep_idx_name:str, multi_mass:float):

    input_Data.pep_idx_name = pep_idx_name
    input_Data.pep_idx_folder = pep_idx_folder
    input_Data.multi_mass = multi_mass

def op_create_peptide_or_spectrum_index(frw82frnw, G, T, mul=1):

    num = int((T - G) * mul) + 1
    index_list = [-1 for i in range(num)]
    for i in range(len(frw82frnw)):
        p = frw82frnw[i]
        index = int((p.mass - G) * mul)
        if index >= num:
            index = num - 1
        if index_list[index] == -1:
            index_list[index] = i
    end_val = len(frw82frnw)
    for i in range(num)[::-1]:
        if index_list[i] == -1:
            index_list[i] = end_val
        else:
            end_val = index_list[i]
    return index_list
