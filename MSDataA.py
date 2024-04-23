class CPeptideIDFeature:
    def __init__(self, score, mean_fra_me, std_fra_me, ion_ratio, inten_ratio, continue_b_score, continue_y_score, continue_score):
        self.score = score
        self.mean_fra_me = mean_fra_me
        self.std_fra_me = std_fra_me
        self.ion_ratio = ion_ratio
        self.inten_ratio = inten_ratio
        self.continue_b_score = continue_b_score
        self.continue_y_score = continue_y_score
        self.continue_score = continue_score

    def get_string(self, split_flag = '\t', fid_flag = False, fid_start = 0):
        if fid_flag:
            return str(fid_start) + ":" + "%.6f"%(self.score) + split_flag + str(fid_start+1) + ":" + "%.6f"%(self.mean_fra_me) + split_flag \
            + str(fid_start+2) + ":" + "%.6f"%(self.std_fra_me) + split_flag + str(fid_start+3) + ":" + "%.6f"%(self.ion_ratio) + split_flag \
            + str(fid_start+4) + ":" + "%.6f"%(self.inten_ratio) + split_flag + str(fid_start+5) + ":" + "%.6f"%(self.continue_b_score) + split_flag \
            + str(fid_start+6) + ":" + "%.6f"%(self.continue_y_score) + split_flag + str(fid_start+7) + ":" + "%.6f"%(self.continue_score)
        return "%.6f"%(self.score) + split_flag + "%.6f"%(self.mean_fra_me) + split_flag + "%.6f"%(self.std_fra_me) + split_flag \
        + "%.6f"%(self.ion_ratio) + split_flag + "%.6f"%(self.inten_ratio) + split_flag \
        + "%.6f"%(self.continue_b_score) + split_flag + "%.6f"%(self.continue_y_score) + split_flag \
        + "%.6f"%(self.continue_score)

    def _print(self):
        print("+++++++++print single score+++++++++++++")
        print("Score:", self.score)
        print("Mean fragment mass error:", self.mean_fra_me)
        print("Std fragment mass error:", self.std_fra_me)
        print("Ion ratio", self.ion_ratio)
        print("Inten ratio", self.inten_ratio)
        print("B continue score", self.continue_b_score)
        print("Y continue score", self.continue_y_score)
        print("B\Y continue score", self.continue_score)

class CPeptideIDInfo:
    def __init__(self, sq:str, mass:float, mods:list=[], pro_list:list=[]):
        self.sq:str = sq
        self.mass:float = mass
        self.mods:list = mods
        self.pro_list:list = pro_list

class CLinkPeptideIDInfo:

    def __init__(self, alpha_peptide:CPeptideIDInfo, beta_peptide:CPeptideIDInfo, alpha_site:int, beta_site:int, score:float=0.0, alpha_score:CPeptideIDFeature=None, beta_score:CPeptideIDFeature=None):

        self.alpha_peptide:CPeptideIDInfo = alpha_peptide
        self.beta_peptide:CPeptideIDInfo = beta_peptide
        self.alpha_site:int = alpha_site
        self.beta_site:int = beta_site
        self.score:float = score
        self.alpha_score:CPeptideIDFeature = alpha_score
        self.beta_score:CPeptideIDFeature = beta_score
        self.svm_score:float = 0.0

class CCrosslinkPSMID:

    def __init__(self, title:str, link_peptide:CLinkPeptideIDInfo, decoy_flag:bool, ranker:int, delta_score:float):
        self.title:str = title
        self.link_peptide:CLinkPeptideIDInfo = link_peptide
        self.decoy_flag:bool = decoy_flag
        self.ranker:int = ranker
        self.delta_score:float = delta_score

class CSinglePeptideIDInfo:

    def __init__(self, peptide_sq: str, mass: float, open_mass: float, pro_flag:int, score: float = 0.0):
        self.peptide_sq = peptide_sq
        self.mass = mass
        self.open_mass = open_mass
        self.pro_flag = pro_flag
        self.score = score

class CSinglePSMID:

    def __init__(self, title:str, peptide:CSinglePeptideIDInfo, decoy_flag:bool, ranker:int, delta_score:float):
        self.title:str = title
        self.peptide:CSinglePeptideIDInfo = peptide
        self.decoy_flag:bool = decoy_flag
        self.ranker:int = ranker
        self.delta_score:float = delta_score