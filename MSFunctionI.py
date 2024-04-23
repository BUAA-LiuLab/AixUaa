from MSSystem import NAME_UAA_PROTEIN, ATOM_MASS_P
from MSDataC import CDataPack
from MSDataPeptide import CInfo1, CProtein, CPeptideInfo
from MSOperatorA import op_FillCInfo1
from MSDataD import CSearchFromBeta
from MSDataB import CPSM
from MSOperatorB import op_FillCSearchFromBeta
from MSFunctionF import CFunctionX12, CFunctionX14, CFunctionX16
from MSFunctionG import CFunctionX18, CFunctionX19
from MSFunctionC import CFunctionX11
from MSFunctionB import CFunctionX1

from MSTool import tool_get_start_T

import operator
import os

class CFunctionX25:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP
        self.functionSpectrum = CFunctionX14(self.dp)

    def SearchMGFs(self, beta_proteins:CProtein, beta_peptide_list:list):

        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        for path_mgf_file in ofir382:
            psm_list = self.SearchOneMGF(path_mgf_file, beta_proteins, beta_peptide_list)
            cmpfun = operator.attrgetter('peptide.score')
            psm_list.sort(key=cmpfun, reverse=True)
            if path_mgf_file.find('/') >= 0:
                result_file_name = path_mgf_file[path_mgf_file.rfind('/') + 1:]
            else:
                result_file_name = path_mgf_file[path_mgf_file.rfind('\\') + 1:]
            result_file_name = result_file_name[:result_file_name.find('.')]
            pkl_file_name = result_file_name + ".pkl"
            txt_file_name = result_file_name + ".txt"
            can_txt_file_name = result_file_name + "_cand.txt"
            print("[Info]Results for {0} are saved to {1}.".format(path_mgf_file, txt_file_name))
            # pickleFunc.save_psm_to_pkl(psm_list, os.path.join(conf.result_folder, pkl_file_name))
            CFunctionX16.WritePSMs(psm_list, os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, txt_file_name), self.dp.myCFG)
            CFunctionX16.WritePSMs(psm_list, os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, can_txt_file_name), self.dp.myCFG, True)
            CFunctionX16.WriteOnepLabel(psm_list, path_mgf_file, self.dp.myCFG.A5_PATH_RESULT_EXPORT, self.dp)

    def SearchOneMGF(self, path_mgf_file:str, beta_proteins:CProtein, beta_peptide_list:list):
        spectrum_list = self.functionSpectrum.LoadSpectrum(path_mgf_file, self.dp.myCFG.D18_PATH_PFIND_RESULT)
        cmpfun = operator.attrgetter('mass')
        spectrum_list.sort(key=cmpfun)
        # 一个mgf文件中所有谱图按质量排好序
        print("[Info]#Spectra in {0} is {1}".format(path_mgf_file, len(spectrum_list)))
        spec_index = self.functionSpectrum.GSMI(spectrum_list, 0, spectrum_list[-1].mass)
        final_res = [[] for i in range(len(spectrum_list))]
        for beta_cewur82214849 in beta_peptide_list:
            # beta_peptide_sq = beta_proteins.sq[beta_cewur82214849.pro_index][beta_cewur82214849.start_pos: beta_cewur82214849.end_pos]
            list_alpha_search_info = self.__captfjewoiuhwqoijfoiwqewqewq(beta_proteins, beta_cewur82214849, spectrum_list, spec_index)
            res = CFunctionX18.fh9ewuhr93hr983h98fwehf9ewh983qe(list_alpha_search_info, self.dp.myCFG.D2_TYPE_THREAD, self.dp.myCFG.D1_NUMBER_THREAD)

            assert (len(res) == len(list_alpha_search_info))
            for i, r in enumerate(res):
                ssi = list_alpha_search_info[i].spectrum_part_start_idx
                for j, spec_r in enumerate(r):
                    if len(final_res[ssi + j]) == 0:
                        final_res[ssi + j] = spec_r
                    else:
                        final_res[ssi+j] += spec_r
                        cmpfun = operator.attrgetter('score')
                        final_res[ssi + j].sort(key=cmpfun, reverse=True)
                    # final_res[ssi+j] = final_res[ssi+j][:self.dp.myCFG.D11_NUMBER_TOP_RESULT]

        psm_list = []
        for i in range(len(final_res)):
            p_list = final_res[i]
            if len(p_list) == 0:
                continue
            best_p = p_list[0]
            cur_spec = spectrum_list[i]
            cur_spec.peaks = []
            psm_list.append(CPSM(cur_spec, best_p, p_list))
        cmpfun = operator.attrgetter('peptide.score')
        psm_list.sort(key=cmpfun, reverse=True)
        return psm_list

    def __captfjewoiuhwqoijfoiwqewqewq(self, beta_proteins:CProtein, beta_cewur82214849:CPeptideInfo, spectrum_list:list, spec_index:list):

        multi = self.dp.myCFG.D12_MULTI_MASS
        beta_mass = beta_cewur82214849.mass

        functionX1 = CFunctionX1(self.dp)
        list_pep_index_tuple = functionX1.GetPepIndex()

        list_search_alpha_info = []
        for i in range(len(list_pep_index_tuple)):
            alpha_G, alpha_T = list_pep_index_tuple[i][2], list_pep_index_tuple[i][3]
            spec_G, spec_T = alpha_G + beta_mass + ATOM_MASS_P, alpha_T + beta_mass + ATOM_MASS_P
            if self.dp.myCFG.C1_TYPE_TOL_PRECURSOR:
                spec_G -= (spec_G * self.dp.myCFG.C2_PPM_TOL_PRECURSOR)
                spec_T += (spec_T * self.dp.myCFG.C2_PPM_TOL_PRECURSOR)
            else:
                spec_G -= self.dp.myCFG.C2_PPM_TOL_PRECURSOR
                spec_T += self.dp.myCFG.C2_PPM_TOL_PRECURSOR
            if spec_G > spectrum_list[-1].mass:
                continue
            if spec_G < 0:
                spec_G = 0.0
            if spec_T > spectrum_list[-1].mass:
                spec_T = spectrum_list[-1].mass
            s_start_ind = spec_index[int(spec_G * multi)]
            while s_start_ind < len(spectrum_list) and spectrum_list[s_start_ind].mass < spec_G:
                s_start_ind += 1
            s_end_ind = spec_index[int(spec_T * multi)]
            while s_end_ind < len(spectrum_list) and spectrum_list[s_end_ind].mass <= spec_T:
                s_end_ind += 1
            one_spectrum_list = spectrum_list[s_start_ind:s_end_ind]
            if len(one_spectrum_list) == 0:
                continue

            search_alpha_info = CSearchFromBeta()
            op_FillCSearchFromBeta(search_alpha_info, one_spectrum_list, beta_proteins, beta_cewur82214849, self.dp, list_pep_index_tuple[i][0], list_pep_index_tuple[i][1], s_start_ind)
            list_search_alpha_info.append(search_alpha_info)

        return list_search_alpha_info


    def __captafmoewiWIEOJRQI3OE2(self, spectrum_list, res, list_alpha_search_info):

        final_res = [[] for i in range(len(spectrum_list))]
        assert (len(res) == len(list_alpha_search_info))
        for i, r in enumerate(res):
            ssi = list_alpha_search_info[i][1]
            for j, spec_r in enumerate(r):
                if len(final_res[ssi + j]) == 0:
                    final_res[ssi + j] = spec_r
                else:
                    final_res[ssi + j] += spec_r
                    cmpfun = operator.attrgetter('score')
                    final_res[ssi + j].sort(key=cmpfun, reverse=True)
                    # final_res[ssi+j] = final_res[ssi+j][:self.dp.myCFG.D11_NUMBER_TOP_RESULT]
        return final_res


class CFunctionX25CLV:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP
        self.functionSpectrum = CFunctionX14(self.dp)

    def SearchMGFs(self, beta_proteins:CProtein, beta_peptide_list:list):

        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        for path_mgf_file in ofir382:
            psm_list = self.SearchOneMGF(path_mgf_file, beta_proteins, beta_peptide_list)
            cmpfun = operator.attrgetter('peptide.alpha_score.score')
            psm_list.sort(key=cmpfun, reverse=True)
            if path_mgf_file.find('/') >= 0:
                result_file_name = path_mgf_file[path_mgf_file.rfind('/') + 1:]
            else:
                result_file_name = path_mgf_file[path_mgf_file.rfind('\\') + 1:]
            result_file_name = result_file_name[:result_file_name.find('.')]
            pkl_file_name = result_file_name + ".pkl"
            txt_file_name = result_file_name + ".txt"
            can_txt_file_name = result_file_name + "_cand.txt"
            print("[Info]Results for {0} are saved to {1}.".format(path_mgf_file, txt_file_name))
            # pickleFunc.save_psm_to_pkl(psm_list, os.path.join(conf.result_folder, pkl_file_name))
            CFunctionX16.WritePSMs(psm_list, os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, txt_file_name), self.dp.myCFG)
            CFunctionX16.WritePSMs(psm_list, os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, can_txt_file_name), self.dp.myCFG, True)
            CFunctionX16.WriteOnepLabel(psm_list, path_mgf_file, self.dp.myCFG.A5_PATH_RESULT_EXPORT, self.dp)

    def SearchOneMGF(self, path_mgf_file:str, beta_proteins:CProtein, beta_peptide_list:list):
        spectrum_list = self.functionSpectrum.LoadSpectrum(path_mgf_file, self.dp.myCFG.D18_PATH_PFIND_RESULT)
        cmpfun = operator.attrgetter('mass')
        spectrum_list.sort(key=cmpfun)
        print("[Info]#Spectra in {0} is {1}".format(path_mgf_file, len(spectrum_list)))
        spec_index = self.functionSpectrum.GSMI(spectrum_list, 0, spectrum_list[-1].mass)

        final_res = [[] for i in range(len(spectrum_list))]
        for beta_cewur82214849 in beta_peptide_list:
            # beta_peptide_sq = beta_proteins.sq[beta_cewur82214849.pro_index][beta_cewur82214849.start_pos: beta_cewur82214849.end_pos]
            list_alpha_search_info = self.__captainDQJWJDWQOIJOISOJIJCSA(beta_proteins, beta_cewur82214849, spectrum_list, spec_index)
            res = CFunctionX19.jfewijfi93wqhjr9i3jeirer(list_alpha_search_info, self.dp.myCFG.D2_TYPE_THREAD, self.dp.myCFG.D1_NUMBER_THREAD)

            assert (len(res) == len(list_alpha_search_info))
            for i, r in enumerate(res):
                ssi = list_alpha_search_info[i].spectrum_part_start_idx
                for j, spec_r in enumerate(r):
                    if len(final_res[ssi + j]) == 0:
                        final_res[ssi + j] = spec_r
                    else:
                        final_res[ssi + j] += spec_r
                    # cmpfun = operator.attrgetter('match_score')
                    # final_res[ssi + j].sort(key=cmpfun, reverse=True)
                    final_res[ssi + j].sort(key=lambda x:(x.match_score, x.TD, -(x.abs_delta_da)), reverse=True)
                    final_res[ssi+j] = final_res[ssi+j][:self.dp.myCFG.D11_NUMBER_TOP_RESULT]

        psm_list = []
        for i in range(len(final_res)):
            p_list = final_res[i]
            if len(p_list) == 0:
                continue
            best_p = p_list[0]

            cur_spec = spectrum_list[i]
            cur_spec.peaks = []
            psm_list.append(CPSM(cur_spec, best_p, p_list))
        cmpfun = operator.attrgetter('peptide.match_score')
        psm_list.sort(key=cmpfun, reverse=True)
        return psm_list

    def __captainDQJWJDWQOIJOISOJIJCSA(self, beta_proteins:CProtein, beta_cewur82214849:CPeptideInfo, spectrum_list:list, spec_index:list):

        multi = self.dp.myCFG.D12_MULTI_MASS
        beta_mass = beta_cewur82214849.mass
        # beta肽质量

        functionX1 = CFunctionX1(self.dp)
        list_pep_index_tuple = functionX1.GetPepIndex()

        list_search_alpha_info = []
        for i in range(len(list_pep_index_tuple)):
            alpha_G, alpha_T = list_pep_index_tuple[i][2], list_pep_index_tuple[i][3]
            spec_G, spec_T = alpha_G + beta_mass + ATOM_MASS_P, alpha_T + beta_mass + ATOM_MASS_P
            if self.dp.myCFG.C1_TYPE_TOL_PRECURSOR:
                spec_G -= (spec_G * self.dp.myCFG.C2_PPM_TOL_PRECURSOR)
                spec_T += (spec_T * self.dp.myCFG.C2_PPM_TOL_PRECURSOR)
            else:
                spec_G -= self.dp.myCFG.C2_PPM_TOL_PRECURSOR
                spec_T += self.dp.myCFG.C2_PPM_TOL_PRECURSOR
            if spec_G > spectrum_list[-1].mass:
                continue
            if spec_G < 0:
                spec_G = 0.0
            if spec_T > spectrum_list[-1].mass:
                spec_T = spectrum_list[-1].mass
            s_start_ind = spec_index[int(spec_G * multi)]
            while s_start_ind < len(spectrum_list) and spectrum_list[s_start_ind].mass < spec_G:
                s_start_ind += 1
            s_end_ind = spec_index[int(spec_T * multi)]
            while s_end_ind < len(spectrum_list) and spectrum_list[s_end_ind].mass <= spec_T:
                s_end_ind += 1
            one_spectrum_list = spectrum_list[s_start_ind:s_end_ind]
            if len(one_spectrum_list) == 0:
                continue

            search_alpha_info = CSearchFromBeta()
            op_FillCSearchFromBeta(search_alpha_info, one_spectrum_list, beta_proteins, beta_cewur82214849, self.dp, list_pep_index_tuple[i][0], list_pep_index_tuple[i][1], s_start_ind)
            list_search_alpha_info.append(search_alpha_info)

        return list_search_alpha_info


    def __caFJQOIWDOIWQSJZLKDt(self, spectrum_list, res, list_alpha_search_info):

        final_res = [[] for i in range(len(spectrum_list))]
        assert (len(res) == len(list_alpha_search_info))
        for i, r in enumerate(res):
            ssi = list_alpha_search_info[i][1]
            for j, spec_r in enumerate(r):
                if len(final_res[ssi + j]) == 0:
                    final_res[ssi + j] = spec_r
                else:
                    final_res[ssi + j] += spec_r
                    cmpfun = operator.attrgetter('score')
                    final_res[ssi + j].sort(key=cmpfun, reverse=True)
        return final_res