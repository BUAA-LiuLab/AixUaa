from MSDataPeptide import CProtein, CPeptideInfo
from MSDataC import CDataPack
from MSDataD import CSearchFromBeta
from MSDataB import CPSM
from MSOperatorB import op_FillCSearchFromBeta
from MSSystem import ATOM_MASS_P

from MSFunctionG import CFunctionX20

from MSFunctionF import CFunctionX14, CFunctionX16

import operator
import os

class CFunctionX27:

    def __init__(self, inputDP:CDataPack):

        self.dp = inputDP
        self.functionSpectrum = CFunctionX14(self.dp)

    def SearchMGFs(self, proteins:CProtein, peptide_list):

        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        for path_mgf_file in ofir382:
            psm_list = self.SearchOneMGF(path_mgf_file, proteins, peptide_list)
            cmpfun = operator.attrgetter('peptide.score')
            psm_list.sort(key=cmpfun, reverse=True)
            if path_mgf_file.find('/') >= 0:
                result_file_name = path_mgf_file[path_mgf_file.rfind('/') + 1:]
            else:
                result_file_name = path_mgf_file[path_mgf_file.rfind('\\') + 1:]
            result_file_name = result_file_name[:result_file_name.find('.')]
            txt_file_name = result_file_name + "_singlePeptide.txt"
            print("[Info]Results for {0} are saved to {1}.".format(path_mgf_file, txt_file_name))
            CFunctionX16.WriteSinglePSMs(psm_list, os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, txt_file_name))  # д�����ļ���

    def SearchOneMGF(self, path_mgf_file:str, beta_proteins:CProtein, beta_peptide_list:list):

        spectrum_list = self.functionSpectrum.LoadSpectrum(path_mgf_file, self.dp.myCFG.D18_PATH_PFIND_RESULT)
        cmpfun = operator.attrgetter('mass')
        spectrum_list.sort(key=cmpfun)
        print("[Info]#Spectra in {0} is {1}".format(path_mgf_file, len(spectrum_list)))
        spec_index = self.functionSpectrum.GSMI(spectrum_list, 0, spectrum_list[-1].mass)
        final_res = [[] for i in range(len(spectrum_list))]
        list_open_search_beta_info = []
        for beta_cewur82214849 in beta_peptide_list:
            open_search_beta_info = self.__captainDMQWOIWQOSASJDOAWJDWQOIJFWQ(beta_proteins, beta_cewur82214849, spectrum_list, spec_index)
            list_open_search_beta_info.append(open_search_beta_info)
        res = CFunctionX20.qwe2qesadfsdfdzgwgewq(list_open_search_beta_info, self.dp.myCFG.D2_TYPE_THREAD, self.dp.myCFG.D1_NUMBER_THREAD)
        assert (len(res) == len(list_open_search_beta_info))
        for i, r in enumerate(res):
            ssi = list_open_search_beta_info[i].spectrum_part_start_idx
            for j, spec_r in enumerate(r):
                final_res[ssi + j].append(spec_r)
                cmpfun = operator.attrgetter('score')
                final_res[ssi + j].sort(key=cmpfun, reverse=True)

        psm_list = []
        for i in range(len(final_res)):
            p_list = final_res[i]
            if len(p_list) == 0: continue
            best_p = p_list[0]
            cur_spec = spectrum_list[i]
            cur_spec.peaks = []
            psm_list.append(CPSM(cur_spec, best_p, p_list))
        cmpfun = operator.attrgetter('peptide.score')
        psm_list.sort(key=cmpfun, reverse=True)

        return psm_list

    def __captainDMQWOIWQOSASJDOAWJDWQOIJFWQ(self, beta_proteins:CProtein, beta_cewur82214849:CPeptideInfo, spectrum_list:list, spec_index:list):

        beta_mass = beta_cewur82214849.mass
        spec_G, spec_T = beta_mass + ATOM_MASS_P - self.dp.myCFG.D17_MASS_WINDOW_BETA, beta_mass + ATOM_MASS_P + self.dp.myCFG.D17_MASS_WINDOW_BETA  # beta�������������Χ
        if spec_G > spectrum_list[-1].mass: return [], 0, -1
        if spec_G < 0: spec_G = 0.0
        if spec_T > spectrum_list[-1].mass: spec_T = spectrum_list[-1].mass

        multi = self.dp.myCFG.D12_MULTI_MASS
        s_start_ind = spec_index[int(spec_G * multi)]
        while s_start_ind < len(spectrum_list) and spectrum_list[s_start_ind].mass < spec_G:
            s_start_ind += 1
        s_end_ind = spec_index[int(spec_T * multi)]
        while s_end_ind < len(spectrum_list) and spectrum_list[s_end_ind].mass <= spec_T:
            s_end_ind += 1


        one_spectrum_list = spectrum_list[s_start_ind:s_end_ind]
        open_search_beta_info = CSearchFromBeta()
        op_FillCSearchFromBeta(open_search_beta_info, one_spectrum_list, beta_proteins, beta_cewur82214849, self.dp, None, None, s_start_ind)

        return open_search_beta_info
