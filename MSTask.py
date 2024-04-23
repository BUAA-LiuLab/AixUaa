from MSDataC import CDataPack
from MSDataPeptide import CInfo1, CProtein

from MSFunctionA import CFunctionConfigIO
from MSFunctionE import CFunctionX4

from MSFunctionF import CFunctionX12, CFunctionX13
from MSSystem import NAME_PROTEIN_PKL, INPUT_SVM_FOLDER
from MSSystem import NAME_UAA_PROTEIN
from MSFunctionC import CFunctionX10, CFunctionX11
from MSOperatorA import op_FillCInfo1
from MSFunctionD import CFunctionX5, CFunctionX6
from MSFunctionI import CFunctionX25, CFunctionX25CLV
from MSFunctionH import CFunctionX23

import os


class CTaskA:
    #建立config文件
    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self, path):
        CFunctionConfigIO.config2file(path, self.dp.myCFG)


class CTaskG:
    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP
    def work(self, argv):
        CFunctionConfigIO.file2config(argv[1], self.dp.myCFG)


class CTaskB:
    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP
    def work(self):
        functionX4 = CFunctionX4(self.dp)
        functionX4.do()


class CTaskH:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def work(self, print_flag=False):

        einfo = CInfo1()
        op_FillCInfo1(einfo, self.dp.myCFG.B1_NAME_ENZYME, self.dp.myCFG.D9_LEN_PEP_UP, self.dp.myCFG.D8_LEN_PEP_LOW, self.dp.myCFG.D7_MASS_PEP_UP, self.dp.myCFG.D6_MASS_PEP_LOW, self.dp.myCFG.B2_TYPE_DIGEST, self.dp.myCFG.B3_NUMBER_MAX_MISS_CLV)

        functionX13 = CFunctionX13()
        functionX10 = CFunctionX10(self.dp, einfo)
        if not (os.path.exists(os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, NAME_PROTEIN_PKL))):
            db_proteins = CFunctionX12.LoadP(self.dp.myCFG.A3_PATH_FASTA)
            functionX10.GMP(db_proteins, self.dp.myCFG.B4_NAME_MOD_FIX, self.dp.myCFG.B5_NAME_MOD_VAR, print_Flag=print_flag)
            functionX13.SavePkl(db_proteins, self.dp.myCFG.A4_PATH_FASTA_EXPORT, NAME_PROTEIN_PKL)


class CTaskI:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def search(self):

        beta_proteins = CProtein()
        list_beta_protein_sq = self.dp.myCFG.B7_UAA_SEQ.split(";")
        for beta_idx, beta_sq in enumerate(list_beta_protein_sq):
            CFunctionX12.APr(beta_proteins, beta_sq, f"{NAME_UAA_PROTEIN}_%d" % beta_idx, "", False)

        beta_peptide_list = self.__captainGetBetaPeptideList(beta_proteins)
        functionSearchForOpenUAALink = CFunctionX25(self.dp)
        functionSearchForOpenUAALink.SearchMGFs(beta_proteins, beta_peptide_list)

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


class CTaskF:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def search(self):
        beta_proteins = CProtein()
        list_beta_protein_sq = self.dp.myCFG.B7_UAA_SEQ.split(";")
        for beta_idx, beta_sq in enumerate(list_beta_protein_sq):
            CFunctionX12.APr(beta_proteins, beta_sq, f"{NAME_UAA_PROTEIN}_%d" % beta_idx, "", False)

        beta_peptide_list = self.__captainGetBetaPeptideList(beta_proteins)
        functionSearchForOpenUAALinkCLV = CFunctionX25CLV(self.dp)
        functionSearchForOpenUAALinkCLV.SearchMGFs(beta_proteins, beta_peptide_list)

    def __captainGetBetaPeptideList(self, beta_proteins: CProtein):

        beta_enzyme_info = CInfo1()
        op_FillCInfo1(beta_enzyme_info, self.dp.myCFG.B1_NAME_ENZYME, self.dp.myCFG.B10_UAA_LEN_UP,
                      self.dp.myCFG.B9_UAA_LEN_LOW, self.dp.myCFG.B12_UAA_MASS_UP, self.dp.myCFG.B11_UAA_MASS_LOW,
                      self.dp.myCFG.B2_TYPE_DIGEST, self.dp.myCFG.B18_UAA_NUMBER_MAX_MISS_CLV)

        functionGenerateBetaPeptide = CFunctionX11(self.dp, beta_enzyme_info)

        beta_peptide_list = functionGenerateBetaPeptide.EPIWJGEWIO(beta_proteins, self.dp.myCFG.B13_UAA_NAME_MOD_FIX, self.dp.myCFG.B14_UAA_NAME_MOD_VAR)
        valid_beta_peptide_list = []
        for beta_cewur82214849 in beta_peptide_list:
            beta_peptide_sq = beta_proteins.sq[beta_cewur82214849.pro_index][beta_cewur82214849.start_pos: beta_cewur82214849.end_pos]
            if beta_peptide_sq.find(self.dp.myCFG.B8_UAA_AA) >= 0:
                valid_beta_peptide_list.append(beta_cewur82214849)
        return valid_beta_peptide_list


class CTaskE:

    def __init__(self, inputDP):
        self.dp = inputDP

    def __captainRerankResult(self):
        functionGetPSMResult = CFunctionX5(self.dp)
        funtionRerank = CFunctionX23(self.dp)
        result_file_list = functionGetPSMResult.FindResFile()
        psm_list = []
        for file_path in result_file_list:
            psm_list += functionGetPSMResult.LoadOpenUAACrosslinkPSM(file_path)
        if not os.path.exists(os.path.join(os.getcwd() + INPUT_SVM_FOLDER)):
            os.makedirs(os.path.join(os.getcwd() + INPUT_SVM_FOLDER))
        feature_write_path = os.path.join(os.getcwd() + INPUT_SVM_FOLDER, functionGetPSMResult.GetTaskId())
        self.dp.SVM_status = funtionRerank.Rerank(psm_list, result_file_list, feature_write_path)

    def RerankResult(self):
        self.__captainRerankResult()


class CTaskD:

    def __init__(self, inputDP):

        self.dp = inputDP
        self.svm_status = self.dp.SVM_status

    def compute_fdr(self):
        functionGetPSMresult = CFunctionX5(self.dp)

        cand_file_list = functionGetPSMresult.FindResFile(rerank_flag=self.svm_status)
        print("[Info] Candidate file list", cand_file_list)
        psm_list = []
        for file_path in cand_file_list:
            psm_list += self._load_result(file_path)
        fdr_score = self._compute_fdr_score(psm_list)
        print("[Info] Score for fdr<={0} is {1}".format(self.dp.myCFG.E1_FDR_PSM, fdr_score))
        self._write_fdr_result(fdr_score)

    def _load_result(self, path):
        functionGetPSMresult = CFunctionX5(self.dp)
        if self.svm_status:
            psm_cand_list = functionGetPSMresult.LoadOpenUAARerankCrosslinkPSM(path)
        else:
            psm_cand_list = functionGetPSMresult.LoadOpenUAACrosslinkPSM(path)
        psm_list = []
        title_dict = {}
        for one_psm in psm_cand_list:
            if one_psm.title not in title_dict: title_dict[one_psm.title] = []
            title_dict[one_psm.title].append((one_psm.link_peptide, one_psm.decoy_flag, one_psm.ranker))
        for title in title_dict:
            for (pep, decoy, ranker) in title_dict[title]:
                if ranker == 1:
                    if not decoy:
                        psm_list.append((title, pep, decoy, ranker))
                    else:
                        cand_list = title_dict[title]
                        cur_decoy = True
                        for (pep2, decoy2, ranker2) in cand_list:
                            if pep2.score == pep.score and not decoy2:
                                cur_decoy = False
                                break
                        psm_list.append((title, pep, cur_decoy, ranker))
        return psm_list

    def _compute_fdr_score(self, psm_list):
        functionCalculateCrosslinkFDR = CFunctionX6(self.dp)
        if self.svm_status:
            score_t = functionCalculateCrosslinkFDR.calculateFDRvalueBySVM(psm_list)
        else:
            score_t = functionCalculateCrosslinkFDR.calculateFDRvalueByORI(psm_list)
        return score_t

    def _write_fdr_result(self, fdr_score):
        functionGetPSMresult = CFunctionX5(self.dp)
        functionCalculateCrosslinkFDR = CFunctionX6(self.dp)
        taskid = functionGetPSMresult.GetTaskId(need_time=False)
        if self.svm_status:
            res_file_list = functionGetPSMresult.FindResFile(rerank_flag=True)
        else:
            res_file_list = functionGetPSMresult.FindResFile()
        res_plabel_file_list = functionGetPSMresult.FindResFile(label_flag=True)
        functionCalculateCrosslinkFDR.write_FDR_result(res_file_list, res_plabel_file_list, fdr_score, taskid, self.svm_status)