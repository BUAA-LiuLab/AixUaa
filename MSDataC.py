
class CConfig:
    # [data]
    A1_TYPE_MS2 = "mgf"
    A2_PATH_MS2 = " #VIP"
    A3_PATH_FASTA = " #VIP"
    A4_PATH_FASTA_EXPORT = "./protein_index/ #VIP"
    A5_PATH_RESULT_EXPORT = "./result/ #VIP"

    # [biology]
    B1_NAME_ENZYME = "trypsin KR C #use ';' to set multiple enzymes"
    B2_TYPE_DIGEST = "0 #0 for specific; 1 for semi-specific; 2 for non-specific"
    B3_NUMBER_MAX_MISS_CLV = "3"
    B4_NAME_MOD_FIX = " #VIP, use ';' to set multiply fixed modifications"
    B5_NAME_MOD_VAR = "Oxidation[M] #VIP, use ';' to set multiply variable modifications"
    B6_NUMBER_MAX_MOD = "3 #Maximum of variable modification in one peptide sequence (not consider the fixed modifications)"
    B7_UAA_SEQ = " #VIP"
    B8_UAA_AA = "U"
    B9_UAA_LEN_LOW = "4"
    B10_UAA_LEN_UP = "20"
    B11_UAA_MASS_LOW = "4"
    B12_UAA_MASS_UP = "20"
    B13_UAA_NAME_MOD_FIX = ""
    B14_UAA_NAME_MOD_VAR = ""
    B15_UAA_COM = " #VIP"
    B16_UAA_NAME_ENZYME = "trypsin KR C # use ';' to set multiply enzymes"
    B17_UAA_TYPE_DIGEST = "0 #0 for specific; 1 for semi-specific; 2 for non-specific"
    B18_UAA_NUMBER_MAX_MISS_CLV = "0"
    B19_UAA_LINKED_AA = " #VIP, beta peptide is linked with which amino acids in alpha peptide, it can be multiply amino acids (e.g., ACDEF)"

    # [mass spectrometry]
    C1_TYPE_TOL_PRECURSOR = "ppm"
    C2_PPM_TOL_PRECURSOR = "20"
    C3_TYPE_TOL_FRAGMENT = "ppm"
    C4_PPM_TOL_FRAGMENT = "20"
    C5_ACTIVATION_TYPE = "HCD"

    # [performance]
    D1_NUMBER_THREAD = "8"
    D2_TYPE_THREAD = "0 #0 is for multi-process (high speed but use more memory); 1 is for multi-thread (low speed but use less memory)"
    D3_NUMBER_SELECT_PEAK = "200"
    D4_NUMBER_SPECTRUM = "10000"
    D5_LEN_MAX_PROTEIN = "100000"
    D6_MASS_PEP_LOW = "400"
    D7_MASS_PEP_UP = "10000"
    D8_LEN_PEP_LOW = "6"
    D9_LEN_PEP_UP = "100"
    D10_INDEX_SPLIT_MASS = "100 #create one pkl file for each 100Da ([0, 100], [100, 200], ..., [9900, 10000])"
    D11_NUMBER_TOP_RESULT = "10 #output top-10 peptides for each spectrum"
    D12_MULTI_MASS = "1 #use mass hash (mass*MUTLI_MASS) to retrive peptide, spectrum or peak, this value of creating peptide index and searching mgf must be same."
    D13_TYPE_TASK = "1"
    D14_TYPE_FILTER_BETA = "1 #whether to filter spectrum which has not matched ion when matching with beta peptide (default is 1)"
    D15_NUMBER_PEAK_BETA = "1 #when 'TYPE_FILTER_BETA' is set as 1, then this value is valid, \
         it will filter spectrum where the number of matched beta ions is less than 'NUMBER_PEAK_BETA' (default is 1)"
    D16_OPEN_SEARCH_SINGLE = "0 # open search for single peptide (beta), 0 means it don't support open search while 1 means it supports"
    D17_MASS_WINDOW_BETA = "300 #when open search only for single beta peptide, the mass window size of open search is 300 Da"
    D18_PATH_PFIND_RESULT = " #path_mgf of pfind result file; if not exists, it can be empty"

    # [filter]
    E1_FDR_PSM = "0.05"

    # [ini]
    F1_PATH_INI_ELEMENT = "./ini/element.ini"
    F2_PATH_INI_AA = "./ini/aa.ini"
    F3_PATH_INI_MOD = "./ini/modification.ini"

    # [advance]
    G1_CLEAVABLE_FLAG = "0 # special for cleavable UAA, 0 means it don't support while 1 means it supports"
    G2_BETA_UAA_NEW_COM = " # the element composition of UAA after cleave"
    G3_ALPHA_AA_TYPE = " # specify the linked amino acids"
    G4_ALPHA_AA_MASS_CHANGE = "0 # the mass change of linked amino acids after cleave"

    def __init__(self):

        self.A1_TYPE_MS2 = CConfig.A1_TYPE_MS2
        self.A2_PATH_MS2 = CConfig.A2_PATH_MS2
        self.A3_PATH_FASTA = CConfig.A3_PATH_FASTA
        self.A4_PATH_FASTA_EXPORT = CConfig.A4_PATH_FASTA_EXPORT
        self.A5_PATH_RESULT_EXPORT = CConfig.A5_PATH_RESULT_EXPORT

        self.B1_NAME_ENZYME = CConfig.B1_NAME_ENZYME
        self.B2_TYPE_DIGEST = CConfig.B2_TYPE_DIGEST
        self.B3_NUMBER_MAX_MISS_CLV = CConfig.B3_NUMBER_MAX_MISS_CLV
        self.B4_NAME_MOD_FIX = CConfig.B4_NAME_MOD_FIX
        self.B5_NAME_MOD_VAR = CConfig.B5_NAME_MOD_VAR
        self.B6_NUMBER_MAX_MOD = CConfig.B6_NUMBER_MAX_MOD
        self.B7_UAA_SEQ = CConfig.B7_UAA_SEQ
        self.B8_UAA_AA = CConfig.B8_UAA_AA
        self.B9_UAA_LEN_LOW = CConfig.B9_UAA_LEN_LOW
        self.B10_UAA_LEN_UP = CConfig.B10_UAA_LEN_UP
        self.B11_UAA_MASS_LOW = CConfig.B11_UAA_MASS_LOW
        self.B12_UAA_MASS_UP = CConfig.B12_UAA_MASS_UP
        self.B13_UAA_NAME_MOD_FIX = CConfig.B13_UAA_NAME_MOD_FIX
        self.B14_UAA_NAME_MOD_VAR = CConfig.B14_UAA_NAME_MOD_VAR
        self.B15_UAA_COM = CConfig.B15_UAA_COM
        self.B16_UAA_NAME_ENZYME = CConfig.B16_UAA_NAME_ENZYME
        self.B17_UAA_TYPE_DIGEST = CConfig.B17_UAA_TYPE_DIGEST
        self.B18_UAA_NUMBER_MAX_MISS_CLV = CConfig.B18_UAA_NUMBER_MAX_MISS_CLV
        self.B19_UAA_LINKED_AA = CConfig.B19_UAA_LINKED_AA

        self.C1_TYPE_TOL_PRECURSOR = CConfig.C1_TYPE_TOL_PRECURSOR
        self.C2_PPM_TOL_PRECURSOR = CConfig.C2_PPM_TOL_PRECURSOR
        self.C3_TYPE_TOL_FRAGMENT = CConfig.C3_TYPE_TOL_FRAGMENT
        self.C4_PPM_TOL_FRAGMENT = CConfig.C4_PPM_TOL_FRAGMENT
        self.C5_ACTIVATION_TYPE = CConfig.C5_ACTIVATION_TYPE

        self.D1_NUMBER_THREAD = CConfig.D1_NUMBER_THREAD
        self.D2_TYPE_THREAD = CConfig.D2_TYPE_THREAD
        self.D3_NUMBER_SELECT_PEAK = CConfig.D3_NUMBER_SELECT_PEAK
        self.D4_NUMBER_SPECTRUM = CConfig.D4_NUMBER_SPECTRUM
        self.D5_LEN_MAX_PROTEIN = CConfig.D5_LEN_MAX_PROTEIN
        self.D6_MASS_PEP_LOW = CConfig.D6_MASS_PEP_LOW
        self.D7_MASS_PEP_UP = CConfig.D7_MASS_PEP_UP
        self.D8_LEN_PEP_LOW = CConfig.D8_LEN_PEP_LOW
        self.D9_LEN_PEP_UP = CConfig.D9_LEN_PEP_UP
        self.D10_INDEX_SPLIT_MASS = CConfig.D10_INDEX_SPLIT_MASS
        self.D11_NUMBER_TOP_RESULT = CConfig.D11_NUMBER_TOP_RESULT
        self.D12_MULTI_MASS = CConfig.D12_MULTI_MASS
        self.D13_TYPE_TASK = CConfig.D13_TYPE_TASK
        self.D14_TYPE_FILTER_BETA = CConfig.D14_TYPE_FILTER_BETA
        self.D15_NUMBER_PEAK_BETA = CConfig.D15_NUMBER_PEAK_BETA
        self.D16_OPEN_SEARCH_SINGLE = CConfig.D16_OPEN_SEARCH_SINGLE
        self.D17_MASS_WINDOW_BETA = CConfig.D17_MASS_WINDOW_BETA
        self.D18_PATH_PFIND_RESULT = CConfig.D18_PATH_PFIND_RESULT

        self.E1_FDR_PSM = CConfig.E1_FDR_PSM

        self.F1_PATH_INI_ELEMENT = CConfig.F1_PATH_INI_ELEMENT
        self.F2_PATH_INI_AA = CConfig.F2_PATH_INI_AA
        self.F3_PATH_INI_MOD = CConfig.F3_PATH_INI_MOD

        self.G1_CLEAVABLE_FLAG = CConfig.G1_CLEAVABLE_FLAG
        self.G2_BETA_UAA_NEW_COM = CConfig.G2_BETA_UAA_NEW_COM
        self.G3_ALPHA_AA_TYPE = CConfig.G3_ALPHA_AA_TYPE
        self.G4_ALPHA_AA_MASS_CHANGE = CConfig.G4_ALPHA_AA_MASS_CHANGE


class CINI:
    DIC_ELEMENT_MASS = {}
    DIC_MOD = {}
    DIC_AA = {}
    MATRIX_GDM = []
    DIC_SET_MOD = {}
    CLV_UAA_MASS = 0.0
    UAA_LINKERS_MASS = 0.0

    def __init__(self):

        self.DIC_ELEMENT_MASS = CINI.DIC_ELEMENT_MASS
        self.DIC_MOD = CINI.DIC_MOD
        self.DIC_AA = CINI.DIC_AA
        self.MATRIX_GDM = CINI.MATRIX_GDM
        self.DIC_SET_MOD = CINI.DIC_SET_MOD
        self.CLV_UAA_MASS = CINI.CLV_UAA_MASS
        self.UAA_LINKERS_MASS = CINI.UAA_LINKERS_MASS

class CDataPack:
    myCFG = CConfig()
    myINI = CINI()
    RerankTime:int = 5

    def __init__(self):

        self.myCFG = CDataPack.myCFG
        self.myINI = CDataPack.myINI
        self.RerankTime = CDataPack.RerankTime
        self.SVM_status = True
