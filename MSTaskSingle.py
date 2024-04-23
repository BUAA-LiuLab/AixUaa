from MSFunctionF import CFunctionX14, CFunctionX16, CFunctionX12
from MSFunctionD import CFunctionX5, CFunctionX7
from MSFunctionC import CFunctionX11
from MSFunctionJ import CFunctionX27
from MSOperatorA import op_create_peptide_or_spectrum_index
from MSDataC import CDataPack
from MSDataPeptide import CInfo1, CProtein
from MSOperatorA import op_FillCInfo1
from MSSystem import ATOM_MASS_P, NAME_UAA_PROTEIN
from MSFunctionG import CFunctionX20

import os
import operator

class CTaskC:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def search(self):

        beta_proteins = CProtein()
        list_beta_protein_sq = self.dp.myCFG.B7_UAA_SEQ.split(";")
        for beta_idx, beta_sq in enumerate(list_beta_protein_sq):
            CFunctionX12.APr(beta_proteins, beta_sq, f"{NAME_UAA_PROTEIN}_%d" % beta_idx, "", False)

        beta_peptide_list = self.__captainGetBetaPeptideList(beta_proteins)
        functionSearchForOpenBeta = CFunctionX27(self.dp)
        functionSearchForOpenBeta.SearchMGFs(beta_proteins, beta_peptide_list)

    def __captainGetBetaPeptideList(self, beta_proteins:CProtein):

        beta_enzyme_info = CInfo1()
        op_FillCInfo1(beta_enzyme_info, self.dp.myCFG.B1_NAME_ENZYME, self.dp.myCFG.B10_UAA_LEN_UP, self.dp.myCFG.B9_UAA_LEN_LOW, self.dp.myCFG.B12_UAA_MASS_UP, self.dp.myCFG.B11_UAA_MASS_LOW, self.dp.myCFG.B2_TYPE_DIGEST, self.dp.myCFG.B18_UAA_NUMBER_MAX_MISS_CLV)

        functionGenerateBetaPeptide = CFunctionX11(self.dp, beta_enzyme_info)
        beta_peptide_list = functionGenerateBetaPeptide.EPIWJGEWIO(beta_proteins, self.dp.myCFG.B4_NAME_MOD_FIX, self.dp.myCFG.B5_NAME_MOD_VAR)
        valid_beta_peptide_list = []
        for beta_cewur82214849 in beta_peptide_list:
            beta_peptide_sq = beta_proteins.sq[beta_cewur82214849.pro_index][beta_cewur82214849.start_pos: beta_cewur82214849.end_pos]
            if beta_peptide_sq.find(self.dp.myCFG.B8_UAA_AA) >= 0:
                valid_beta_peptide_list.append(beta_cewur82214849)
        return valid_beta_peptide_list


class CTaskJ:

    def __init__(self, inputDP):
        self.dp = inputDP

    def _load_result(self, file_path):
        functionGetPSMResult = CFunctionX5(self.dp)
        psm_list = functionGetPSMResult.LoadOpenUAASinglePSM(file_path)
        return psm_list

    def __captainComputeFDRScore(self, psm_list):
        functionCalculateSingleFDR = CFunctionX7(self.dp)
        psm_list.sort(key=lambda k: k.peptide.score, reverse=True)
        score_t = functionCalculateSingleFDR.CalculateFDR(psm_list)
        return score_t

    def _write_fdr_result(self, fdr_score):
        functionGetPSMResult = CFunctionX5(self.dp)
        functionCalculateSingleFDR = CFunctionX7(self.dp)
        taskid = functionGetPSMResult.GetTaskId(False)
        res_file_list = functionGetPSMResult.FindSingleResFile()
        functionCalculateSingleFDR.write_FDR_result(res_file_list, fdr_score, taskid)

    def compute_fdr(self):
        psm_list = []
        functionGetPSMResult = CFunctionX5(self.dp)
        cand_file_list = functionGetPSMResult.FindSingleResFile()
        print("[Info] Candidate file list", cand_file_list)
        for file_path in cand_file_list:
            psm_list += self._load_result(file_path)
        fdr_score = self.__captainComputeFDRScore(psm_list)
        print("[Info] Score for fdr<={0} is {1}".format(self.dp.myCFG.E1_FDR_PSM, fdr_score))
        self._write_fdr_result(fdr_score)