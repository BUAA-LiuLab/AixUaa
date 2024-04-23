from MSDataPeptide import CPeptideInfo

class CSingleScoreFeature:
    def __init__(self, match_score:float, mean_fra_me:float, std_fra_me:float, ion_ratio:float, inten_ratio:float, continue_b_score:float, continue_y_score:float, continue_score:float, match_peak_info:list, uaa_clv=False):
        self.score:float = match_score + continue_b_score + continue_y_score + continue_score
        self.match_score: float = match_score
        self.mean_fra_me:float = mean_fra_me
        self.std_fra_me:float = std_fra_me
        self.ion_ratio:float = ion_ratio
        self.inten_ratio:float = inten_ratio
        self.continue_b_score:float = continue_b_score
        self.continue_y_score:float = continue_y_score
        self.continue_score:float = continue_score
        self.match_peak_info:list = match_peak_info
        self.uaa_clv:bool = uaa_clv

    def get_string(self, split_flag = '\t', fid_flag = False, fid_start = 0):
        if fid_flag:
            return str(fid_start) + ":" + "%.6f"%(self.match_score) + split_flag + str(fid_start+1) + ":" + "%.6f"%(self.mean_fra_me) + split_flag \
            + str(fid_start+2) + ":" + "%.6f"%(self.std_fra_me) + split_flag + str(fid_start+3) + ":" + "%.6f"%(self.ion_ratio) + split_flag \
            + str(fid_start+4) + ":" + "%.6f"%(self.inten_ratio) + split_flag + str(fid_start+5) + ":" + "%.6f"%(self.continue_b_score) + split_flag \
            + str(fid_start+6) + ":" + "%.6f"%(self.continue_y_score) + split_flag + str(fid_start+7) + ":" + "%.6f"%(self.continue_score)
        return "%.6f"%(self.match_score) + split_flag + "%.6f"%(self.mean_fra_me) + split_flag + "%.6f"%(self.std_fra_me) + split_flag \
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

class CLinkPeptide:
    def __init__(self, alpha_peptide:CPeptideInfo, alpha_site, beta_peptide:CPeptideInfo, beta_site, score:float=0.0, alpha_score:CSingleScoreFeature=None, beta_score:CSingleScoreFeature=None):
        #alpha_score and beta_score are class objects (see "CSingle_Score_Feature")
        self.alpha_peptide:CPeptideInfo = alpha_peptide
        self.beta_peptide:CPeptideInfo = beta_peptide
        self.alpha_site = alpha_site
        self.beta_site = beta_site
        self.score:float = score
        self.match_score:float = 0.0
        self.alpha_score:CSingleScoreFeature = alpha_score
        self.beta_score:CSingleScoreFeature = beta_score
        self.alpha_sq:str = ""
        self.beta_sq:str = ""
        self.TD:int = 1
        self.abs_delta_da:float = 0.0
        self.clv:int = 0

class CPSM:
    def __init__(self, spec, peptide:CLinkPeptide, peptide_list, decoy_flag=False):
        self.spec = spec
        self.peptide = peptide
        self.decoy_flag = decoy_flag
        self.peptide_list = peptide_list

class CSinglePeptide:
    def __init__(self, peptide:CPeptideInfo, peptide_sq:str, open_site:int, open_mass:float, pro_index:int, score:float=0.0):
        self.peptide:CPeptideInfo = peptide
        self.open_site:int = open_site
        self.open_mass:float = open_mass
        self.pro_index:int = pro_index
        self.score:float = score
        self.peptide_sq:str = peptide_sq
