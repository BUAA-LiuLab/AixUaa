MOLECULE_MASS_H2O = 18.0105647
AA_NUM = 26

ATOM_MASS_P = 1.00727647012
AVRG_AA_MASS = 110.0

VALUE_MAX_PROTEIN_NUM = 50000
VALUE_MAX_SCAN = 2000000

NAME_MOD_PKL = "modification.pkl"
NAME_PEPTIDE_PKL = "peptide.pkl"
NAME_PROTEIN_PKL = "protein.pkl"
NAME_UAA_PROTEIN = "UAAProtein"

INPUT_SVM_FOLDER = "\\svm_tool\\search_alpha_info\\"
MODEL_PATH = "\\svm_tool\\models\\v1"

# INPUT_SVM_FOLDER = "./svm_tool/search_alpha_info/"
# MODEL_PATH = "./svm_tool/models/v1"

SINGLE_PSM = "res_single_psm"
SINGLE_PEP = "res_single"
XLINK_PSM = "res_psm"
XLINK_PEP = "res"
SINGLE_RES_PATH = ""

EXPIRATION_TIME = {'Year': 2024, 'Month': 5, 'Day': 30}

SOFTWARE_NAME = "AixUaa"

CONFIG_NAME = "config_example.txt"

INFO_TO_USER_Staff = (
    f'\n[{SOFTWARE_NAME}] {SOFTWARE_NAME} is expired! Please send e-mail to cliu126@126.com for the new version.',
    f'\n[{SOFTWARE_NAME}] Warning! The current license will expired in 7 days. Please send e-mail to cliu126@126.com for the new version.',
    f'\n[{SOFTWARE_NAME}] Starting...',
    f'\n[{SOFTWARE_NAME}] Finished!',
    f'\n[{SOFTWARE_NAME}] Writing config file...',)

INFO_TO_USER_TaskRead = (
    f'\n[{SOFTWARE_NAME}] Reading MS1 file...',
    f'\n[{SOFTWARE_NAME}] Reading MS2 file...',
    f'\n[{SOFTWARE_NAME}] Reading identification results: ',)


INFO_TO_USER_TaskDraw = (
    f'\n[{SOFTWARE_NAME}] Drawing figures...',)

INFO_TO_USER_TaskExport = (
    f'\n[{SOFTWARE_NAME}] Exporting...',)

INFO_TO_USER_TaskXtract = (

    f'\n[{SOFTWARE_NAME}] Xtracting...',
    f'\n[{SOFTWARE_NAME}] ms1 and ms2 files are existed...',)


INFO_TO_USER_TaskProfile = (

    f'\n[{SOFTWARE_NAME}] Creating index file for ms1...',
    f'\n[{SOFTWARE_NAME}] Reconstructing chromatograms...',
    f'\n[{SOFTWARE_NAME}] #Chromatograms: ',)


FILENAME_EXPORT = (
    'INFO_MS1.txt',
    'INFO_MS2.txt',
    'INFO_Cycle.txt',
    'INFO_ID.txt',
    'INFO_Mass_Deviation.txt',
    'INFO_Chromatography.txt',)

FILENAME_DRAW = (
    'Fig_Histogram_IonInjectionTime_MS1.pdf',
    'Fig_Histogram_IonInjectionTime_MS2.pdf',
    'Fig_Histogram_ElutingTime.pdf',)


NAME_FOR_MSREGISTER = ["MSRegisterExportDLL.dll",
                       f"{SOFTWARE_NAME}.license",
                       f"{SOFTWARE_NAME}.probation"]

LINK_APPLICATION = "https://n5wjnec5xs7pr1sh.mikecrm.com/iuEDadM"

NAME_PUBLIC_SECURITY = "MSSource.dll"

