from MSDataC import CDataPack
from MSDataD import CSearchFromBeta, CSpectrum
from MSDataPeptide import CPeptideInfo, CProtein

import copy


def op_FillCSearchFromBeta(input_Data:CSearchFromBeta, spectrum_info_list:list, beta_proteins:CProtein, beta_info:CPeptideInfo, inputDP:CDataPack, name_pep_index:str, name_pep_index_ind:str, spectrum_part_idx:int):

    input_Data.list_spectrum_info = spectrum_info_list
    input_Data.beta_info = beta_info
    input_Data.beta_sq = beta_proteins.sq[beta_info.pro_index][beta_info.start_pos: beta_info.end_pos]
    input_Data.dp = inputDP
    input_Data.name_pep_index = name_pep_index
    input_Data.name_pep_index_ind = name_pep_index_ind
    input_Data.spectrum_part_start_idx = spectrum_part_idx