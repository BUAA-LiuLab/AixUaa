# 这里定义的是产生肽段需要的数据结构

class CInfo1:

    def __init__(self):
        self.digest_type:int = 0
        self.length_up:int = 0
        self.length_low:int = 0
        self.mass_up = 0.0
        self.mass_low = 0.0
        self.miss_clv:int = 0
        self.enzyme_flag:list = []
        self.enzyme_aa:list = []

class CMod:
    def __init__(self, name, mass, type, aa):
        self.name = name
        self.mass = mass
        self.type = type
        self.aa = aa

class CProtein:
    def __init__(self):
        self.ac = []
        self.de = []
        self.sq = []

class CPeptideInfo:
    def __init__(self, start_pos, end_pos, pro_index, pro_len):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.pro_index = pro_index
        self.pro_len = pro_len
        self.mass = 0.0
        self.gdm = 0.0
        self.mods = []
        self.pro_index_list = []

class CPeptidePKL:
    def __init__(self, pro_index, start_pos, end_pos, mass, mods=[], gdm=0.0):
        self.pro_index = pro_index
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.mass = mass
        self.mods = mods
        self.gdm = gdm
        self.pro_index_list = []

class CModSite:
    def __init__(self, mod_name, site, mass):
        self.mod_name = mod_name
        self.site = site
        self.mass = mass

class CAddModInfo:
    def __init__(self):
        self.fix_N_aa = {}
        self.fix_C_aa = {}
        self.fix_aa = {}
        self.fix_PN_aa = {}
        self.fix_PC_aa = {}
        self.var_N_aa = {}
        self.var_C_aa = {}
        self.var_aa = {}
        self.var_PN_aa = {}
        self.var_PC_aa = {}

        self.modname2type = {}
        self.modname2index = {}
        self.mod_list = []
class CPeptideIndexSortInfo:

    def __init__(self):
        self.pep_idx_name:str = ""
        self.pep_idx_folder:str = ""
        self.multi_mass:int = 100
