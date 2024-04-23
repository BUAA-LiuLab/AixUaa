from MSDataPeptide import CPeptideInfo
from MSDataC import CDataPack

class CPeak:
    def __init__(self, mz, inten):
        self.mz = mz
        self.inten = inten
class CSpectrum:
    def __init__(self, title, mass, charge, peaks=[]):
        self.title = title
        self.mass = mass
        self.charge = charge
        self.peaks = peaks

class CSearchFromBeta:

    def __init__(self):

        self.name_pep_index:str= ""
        self.name_pep_index_ind:str= ""
        self.list_spectrum_info:list=None
        self.beta_info:CPeptideInfo=None
        self.beta_sq:str=""
        self.dp:CDataPack=None
        self.spectrum_part_start_idx:int=0