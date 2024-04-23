from MSSystem import NAME_MOD_PKL, ATOM_MASS_P, NAME_PROTEIN_PKL
from MSDataPeptide import CPeptideIndexSortInfo, CProtein, CInfo1
from MSOperatorA import op_FillCInfo1
from MSDataD import CSearchFromBeta, CSpectrum, CPeptideInfo
from MSDataB import CLinkPeptide, CSingleScoreFeature, CSinglePeptide
from MSDataC import CDataPack
from MSTool import toolCheckLinkSiteValid, tool_get_start_T
from MSFunctionF import CFunctionX13
from MSFunctionE import CFunctionX2

from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
import operator
import os
import time

try:
    import cPickle as pickle
except:
    import pickle


class CFunctionX17:

    @staticmethod
    def jdOSFEWODSFOJNFOWFOfwoenfewofew(list_sort_pep_idx_info, type_thread:int=0, number_thread:int=1):
        if type_thread == 0:
            pool = Pool(processes=number_thread)
        else:
            pool = ThreadPool(number_thread, list_sort_pep_idx_info)
        pool.map(CFunctionX17.fjsoifjoiewjro3, list_sort_pep_idx_info)
        pool.close()
        pool.join()

    @staticmethod
    def fjsoifjoiewjro3(sort_pep_idx_info:CPeptideIndexSortInfo):
        G, T = 0.0, 0.0
        start_end_str = sort_pep_idx_info.pep_idx_name
        start_end_str = sort_pep_idx_info.pep_idx_name[:start_end_str.find('.')]
        # 不要后缀
        index_str = start_end_str.find('_')
        # _这个是质量分割符
        # 表示这个pkl文件对应的质量区间
        G = int(start_end_str[:index_str]) * 1.0
        T = int(start_end_str[index_str + 1:]) * 1.0

        frw82frnw = []
        path_pep_idx = os.path.join(sort_pep_idx_info.pep_idx_folder, sort_pep_idx_info.pep_idx_name)
        f = open(path_pep_idx, 'rb')
        while True:
            try:
                one_peptide = pickle.load(f)
                frw82frnw.append(one_peptide)
            except Exception:
                break
        f.close()
        cmpfun = operator.attrgetter('mass', 'gdm')
        frw82frnw.sort(key=cmpfun)
        old_num = len(frw82frnw)
        iuoytt = CFunctionX17.__rwqqvnhf(frw82frnw)
        print("[Info]#Peptides in {0} file: {1} to {2}".format(path_pep_idx, old_num, len(iuoytt)))
        index_list = CFunctionX17.__cfjewiurh3u2h94328u8few9(iuoytt, G, T, sort_pep_idx_info.multi_mass)
        fw = open(path_pep_idx, 'wb')
        pickle.dump(iuoytt, fw)
        fw.close()
        index_pkl_file = path_pep_idx[:path_pep_idx.rfind('.')] + "_ind.pkl"
        fw = open(index_pkl_file, 'wb')
        pickle.dump(index_list, fw)
        fw.close()

    @staticmethod
    def __rwqqvnhf(frw82frnw):
        start_i = 0
        iuoytt = []
        while start_i < len(frw82frnw):
            iuoytt.append(frw82frnw[start_i])
            cur_p = frw82frnw[start_i]
            cur_mass = cur_p.mass
            cur_gdm = cur_p.gdm
            start_j = start_i + 1
            cur_pro_index_list = []
            cur_pro_index_list += cur_p.pro_index_list
            while start_j < len(frw82frnw) and frw82frnw[start_j].mass == cur_mass and frw82frnw[start_j].gdm == cur_gdm:
                cur_pro_index_list += frw82frnw[start_j].pro_index_list
                start_j += 1
            iuoytt[-1].pro_index_list = cur_pro_index_list
            start_i = start_j
        return iuoytt

    @staticmethod
    def __cfjewiurh3u2h94328u8few9(frw82frnw, G, T, mul=1):
        num = int((T - G) * mul) + 1
        index_list = [-1 for i in range(num)]
        for i in range(len(frw82frnw)):
            p = frw82frnw[i]
            index = int((p.mass - G) * mul)
            if index >= num: index = num - 1
            if index_list[index] == -1:
                index_list[index] = i
        end_val = len(frw82frnw)
        for i in range(num)[::-1]:
            if index_list[i] == -1:
                index_list[i] = end_val
            else:
                end_val = index_list[i]
        return index_list


class CFunctionX18:

    @staticmethod
    def fh9ewuhr93hr983h98fwehf9ewh983qe(input_data_list, type_thread:int=0, number_thread:int=1):
        if len(input_data_list) < number_thread:
            number_thread = len(input_data_list)
        if type_thread == 0:
            pool = Pool(processes=number_thread)
        else:
            pool = ThreadPool(number_thread, input_data_list)
        res = pool.map(CFunctionX18.jdsoFJEWIFJEIWfjwieqjr343, input_data_list)
        return res

    @staticmethod
    def jdsoFJEWIFJEIWfjwieqjr343(search_alpha_info:CSearchFromBeta):
        spectrum_list = search_alpha_info.list_spectrum_info
        pkl_file = search_alpha_info.name_pep_index
        index_pkl_file = search_alpha_info.name_pep_index_ind
        beta_cewur82214849 = search_alpha_info.beta_info
        beta_peptide_sq = search_alpha_info.beta_sq
        beta_peptide_mass = beta_cewur82214849.mass
        dp = search_alpha_info.dp

        alpha_enzyme_info = CInfo1()
        op_FillCInfo1(alpha_enzyme_info, dp.myCFG.B1_NAME_ENZYME, dp.myCFG.D9_LEN_PEP_UP,
                      dp.myCFG.D8_LEN_PEP_LOW, dp.myCFG.D7_MASS_PEP_UP, dp.myCFG.D6_MASS_PEP_LOW,
                      dp.myCFG.B2_TYPE_DIGEST, dp.myCFG.B3_NUMBER_MAX_MISS_CLV)

        C_str, N_str = "", ""
        for i in range(len(alpha_enzyme_info.enzyme_aa)):  # 酶解氨基酸位点C端和N端酶
            if alpha_enzyme_info.enzyme_flag[i] == "N":
                N_str += alpha_enzyme_info.enzyme_aa[i]
            elif alpha_enzyme_info.enzyme_flag[i] == "C":
                C_str += alpha_enzyme_info.enzyme_aa[i]

        print("[Info]#Spectra for {0} is {1}".format(pkl_file, len(spectrum_list)))
        rew = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, search_alpha_info.name_pep_index)
        if len(rew) == 0:
            return [[] for i in range(len(spectrum_list))]
        res = []
        rew_ind = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, search_alpha_info.name_pep_index_ind)
        alpha_proteins = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, NAME_PROTEIN_PKL)

        beta_link_site = beta_peptide_sq.find(dp.myCFG.B8_UAA_AA)
        pep_idx_G, pep_idx_T = tool_get_start_T(pkl_file)

        for spectrum in spectrum_list:
            start_time = time.time()
            list_link_cewur82214849 = []
            beta_score_feature = CFunctionX2.EFHU3H2JFOIWEO3EWR32fdewfwwf(dp, spectrum, beta_cewur82214849, beta_peptide_sq, beta_link_site, dp.myCFG.D12_MULTI_MASS)
            candidate_peptide_list = CFunctionX18.__cafjeoiwnfoij3oi4joi3foi(search_alpha_info, pep_idx_G, spectrum, rew, rew_ind)
            if len(candidate_peptide_list) == 0:
                res.append([])
                continue
            for alpha_cewur82214849 in candidate_peptide_list:
                alpha_peptide_sq = alpha_proteins.sq[alpha_cewur82214849.pro_index][alpha_cewur82214849.start_pos: alpha_cewur82214849.end_pos]
                list_alpha_link_site = CFunctionX18.__cjdoqiwjoiqwhjer9iq32121421(alpha_cewur82214849, alpha_proteins, dp.myCFG.B19_UAA_LINKED_AA, C_str, N_str)
                for alpha_link_site in list_alpha_link_site:
                    link_cewur82214849 = CFunctionX18.oEWIJ3RJFsajfiowejr309843(dp, spectrum, alpha_cewur82214849, alpha_link_site, alpha_peptide_sq, beta_cewur82214849, beta_link_site, beta_peptide_sq, beta_score_feature)
                    list_link_cewur82214849.append(link_cewur82214849)
            end_time = time.time()
            cmpfun = operator.attrgetter('score')
            list_link_cewur82214849.sort(key=cmpfun, reverse=True)
            CFunctionX18.__captfkeoiewio3QIEJDEE(list_link_cewur82214849, alpha_proteins)
            list_final_link_cewur82214849 = list_link_cewur82214849[:dp.myCFG.D11_NUMBER_TOP_RESULT]
            res.append(list_final_link_cewur82214849)
        assert (len(res) == len(spectrum_list))
        return res

    @staticmethod
    def __cafjeoiwnfoij3oi4joi3foi(search_alpha_info:CSearchFromBeta, pep_index_G, spectrum:CSpectrum, alpha_rew:list, alpha_rew_ind):
        start_time = time.time()
        alpha_mass = spectrum.mass - search_alpha_info.beta_info.mass - ATOM_MASS_P
        if search_alpha_info.dp.myCFG.C1_TYPE_TOL_PRECURSOR:
            tol_mass = spectrum.mass * search_alpha_info.dp.myCFG.C2_PPM_TOL_PRECURSOR
        else:
            tol_mass = search_alpha_info.dp.myCFG.C2_PPM_TOL_PRECURSOR
        G = alpha_mass - tol_mass
        T = alpha_mass + tol_mass
        if G > alpha_rew[-1].mass:
            return []
        if G < 0:
            G = 0.0
        if T > alpha_rew[-1].mass:
            T = alpha_rew[-1].mass
        s_start_ind = alpha_rew_ind[int((G - pep_index_G) * search_alpha_info.dp.myCFG.D12_MULTI_MASS)]
        while s_start_ind < len(alpha_rew) and alpha_rew[s_start_ind].mass < G:
            s_start_ind += 1
        s_end_ind = alpha_rew_ind[int((T - pep_index_G) * search_alpha_info.dp.myCFG.D12_MULTI_MASS)]
        while s_end_ind < len(alpha_rew) and alpha_rew[s_end_ind].mass <= T:
            s_end_ind += 1
        if s_end_ind - s_start_ind <= 0:
            return []
        candidate_peptide_list = alpha_rew[s_start_ind: s_end_ind]  # 获取候选肽段区间
        end_time = time.time()
        return candidate_peptide_list

    @staticmethod
    def __cjdoqiwjoiqwhjer9iq32121421(alpha_cewur82214849:CPeptideInfo, alpha_proteins:CProtein, set_link_site:str, C_str:str= "", N_str:str= ""):
        alpha_peptide_sq = alpha_proteins.sq[alpha_cewur82214849.pro_index][alpha_cewur82214849.start_pos:alpha_cewur82214849.end_pos]
        alpha_vnbbmx = []
        mod_site_dict = set()
        for m in alpha_cewur82214849.mods:
            mod_site_dict.add(m.site - 1)
        for s in range(len(alpha_peptide_sq)):
            if alpha_peptide_sq[s] not in set_link_site: continue
            if s in mod_site_dict: continue
            pro_N = False
            pro_C = False
            if s == 0:
                if alpha_cewur82214849.start_pos == 0:
                    pro_N = True
            if s == len(alpha_peptide_sq) - 1:
                if alpha_cewur82214849.end_pos == alpha_cewur82214849.pro_len:
                    pro_C = True
            if not toolCheckLinkSiteValid(C_str, N_str, alpha_peptide_sq, s, pro_N, pro_C): continue
            alpha_vnbbmx.append(s)
        return alpha_vnbbmx

    @staticmethod
    def oEWIJ3RJFsajfiowejr309843(dp:CDataPack, spectrum:CSpectrum, alpha_cewur82214849:CPeptideInfo, alpha_link_site:int, alpha_peptide_sq:str, beta_cewur82214849:CPeptideInfo, beta_link_site:int, beta_peptide_sq:str, beta_score_feature:CSingleScoreFeature, charge_list:list=[1, 2], link_mass:float=None):
        link_cewur82214849 = CLinkPeptide(alpha_cewur82214849, alpha_link_site, beta_cewur82214849, beta_link_site)
        alpha_score_feature = CFunctionX2.EFHU3H2JFOIWEO3EWR32fdewfwwf(dp, spectrum, alpha_cewur82214849, alpha_peptide_sq, alpha_link_site, dp.myCFG.D12_MULTI_MASS, charge_list=charge_list, link_mass=link_mass, alpha_or_beta="α")
        link_cewur82214849.score = alpha_score_feature.score
        link_cewur82214849.alpha_score = alpha_score_feature
        link_cewur82214849.beta_score = beta_score_feature
        link_cewur82214849.alpha_sq = alpha_peptide_sq
        link_cewur82214849.beta_sq = beta_peptide_sq
        return link_cewur82214849

    @staticmethod
    def __captfkeoiewio3QIEJDEE(list_link_cewur82214849:list, alpha_proteins:CProtein):

        if len(list_link_cewur82214849) <= 1:
            return list_link_cewur82214849
        ac = alpha_proteins.ac[list_link_cewur82214849[0].alpha_peptide.pro_index]
        if not ac.startswith("REV_"):
            return list_link_cewur82214849
        if list_link_cewur82214849[0].score - list_link_cewur82214849[1].score < 0.0001:
            return list_link_cewur82214849
        for i in range(1, len(list_link_cewur82214849)):
            if list_link_cewur82214849[i].score < list_link_cewur82214849[0].score:
                break
            new_ac = alpha_proteins.ac[0][list_link_cewur82214849[i].alpha_peptide.pro_index]
            if not new_ac.startswith("REV_"):
                break
            if i < len(list_link_cewur82214849) and list_link_cewur82214849[i].score - list_link_cewur82214849[0].score <= 0.0001:
                tmp = list_link_cewur82214849[0]
                list_link_cewur82214849[0] = list_link_cewur82214849[i]
                list_link_cewur82214849[i] = tmp


class CFunctionX19:

    @staticmethod
    def jfewijfi93wqhjr9i3jeirer(input_data_list, type_thread:int=0, number_thread:int=1):
        if len(input_data_list) < number_thread:
            number_thread = len(input_data_list)
        if number_thread == 0:
            number_thread = 1
        if type_thread == 0:
            pool = Pool(processes=number_thread)  # multiply processes
        else:
            pool = ThreadPool(number_thread)  # multiply threads
        res = pool.map(CFunctionX19.jfewijfoiewhjr93qu984ndsojvklokje, input_data_list)

        pool.close()
        pool.join()
        return res

    @staticmethod
    def addtest(idx:int):
        return idx

    @staticmethod
    def jfewijfoiewhjr93qu984ndsojvklokje(search_alpha_info:CSearchFromBeta):
        spectrum_list = search_alpha_info.list_spectrum_info
        pkl_file = search_alpha_info.name_pep_index
        index_pkl_file = search_alpha_info.name_pep_index_ind
        beta_cewur82214849 = search_alpha_info.beta_info
        beta_peptide_sq = search_alpha_info.beta_sq
        beta_peptide_mass = beta_cewur82214849.mass
        dp = search_alpha_info.dp

        alpha_enzyme_info = CInfo1()
        op_FillCInfo1(alpha_enzyme_info, dp.myCFG.B1_NAME_ENZYME, dp.myCFG.D9_LEN_PEP_UP,
                      dp.myCFG.D8_LEN_PEP_LOW, dp.myCFG.D7_MASS_PEP_UP, dp.myCFG.D6_MASS_PEP_LOW,
                      dp.myCFG.B2_TYPE_DIGEST, dp.myCFG.B3_NUMBER_MAX_MISS_CLV)

        C_str, N_str = "", ""
        for i in range(len(alpha_enzyme_info.enzyme_aa)):
            if alpha_enzyme_info.enzyme_flag[i] == "N":
                N_str += alpha_enzyme_info.enzyme_aa[i]
            elif alpha_enzyme_info.enzyme_flag[i] == "C":
                C_str += alpha_enzyme_info.enzyme_aa[i]

        print("[Info]#Spectra for {0} is {1}".format(pkl_file, len(spectrum_list)))
        rew = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, search_alpha_info.name_pep_index)

        if len(rew) == 0:
            return [[] for i in range(len(spectrum_list))]
        res = []
        rew_ind = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, search_alpha_info.name_pep_index_ind)
        alpha_proteins = CFunctionX13.LoadPkl(dp.myCFG.A4_PATH_FASTA_EXPORT, NAME_PROTEIN_PKL)

        beta_link_site = beta_peptide_sq.find(dp.myCFG.B8_UAA_AA)
        pep_idx_G, pep_idx_T = tool_get_start_T(pkl_file)

        for spectrum in spectrum_list:
            start_time = time.time()
            list_link_cewur82214849 = []
            beta_clv_score_feature = CFunctionX2.EFHIUWHFU32Nfuewinwe2134(dp, spectrum, beta_cewur82214849, beta_peptide_sq, beta_link_site, dp.myCFG.D12_MULTI_MASS)
            charge = min(spectrum.charge, 4)
            charge_list = [i for i in range(1, charge)]
            beta_regular_score_feature = CFunctionX2.EFHU3H2JFOIWEO3EWR32fdewfwwf(dp, spectrum, beta_cewur82214849, beta_peptide_sq, beta_link_site, dp.myCFG.D12_MULTI_MASS, charge_list=charge_list)

            candidate_peptide_list = CFunctionX19.__captainmvsijefoijrqo(search_alpha_info, pep_idx_G, spectrum, rew, rew_ind)
            if len(candidate_peptide_list) == 0:
                res.append([])
                continue
            for alpha_cewur82214849 in candidate_peptide_list:
                alpha_peptide_sq = alpha_proteins.sq[alpha_cewur82214849.pro_index][alpha_cewur82214849.start_pos: alpha_cewur82214849.end_pos]
                list_alpha_link_site = CFunctionX19.__capjfewij0iwji3j340apowqpokdw(alpha_cewur82214849, alpha_proteins, dp.myCFG.B19_UAA_LINKED_AA, C_str, N_str)
                list_alpha_clv_link_site = []
                for alpha_link_site in list_alpha_link_site:
                    if alpha_peptide_sq[alpha_link_site] in dp.myCFG.G3_ALPHA_AA_TYPE:
                        list_alpha_clv_link_site.append(alpha_link_site)
                        continue
                    if len(beta_regular_score_feature.match_peak_info) < dp.myCFG.D15_NUMBER_PEAK_BETA:
                        continue
                    link_cewur82214849 = CFunctionX19.Gefewjiojfoiewjroi3wq(dp, spectrum, alpha_cewur82214849, alpha_link_site, alpha_peptide_sq, beta_cewur82214849, beta_link_site, beta_peptide_sq, beta_regular_score_feature, charge_list=charge_list)
                    if link_cewur82214849.alpha_score.match_score < 1e-6:
                        continue
                    list_link_cewur82214849.append(link_cewur82214849)
                if len(list_alpha_clv_link_site) > 0:
                    if len(beta_clv_score_feature.match_peak_info) < dp.myCFG.D15_NUMBER_PEAK_BETA:
                        continue
                    link_cewur82214849 = CFunctionX19.Gefewjiojfoiewjroi3wq(dp, spectrum, alpha_cewur82214849, list_alpha_clv_link_site[0], alpha_peptide_sq, beta_cewur82214849, beta_link_site, beta_peptide_sq, beta_clv_score_feature)
                    if link_cewur82214849.alpha_score.match_score < 1e-6:
                        continue
                    link_cewur82214849.alpha_site = list_alpha_clv_link_site
                    list_link_cewur82214849.append(link_cewur82214849)
            end_time = time.time()
            # print("[Debug]Time match for spectrum {0} is {1}".format(spectrum.title, (end_time - start_time)))
            cmpfun = operator.attrgetter('match_score')
            list_link_cewur82214849.sort(key=cmpfun, reverse=True)
            CFunctionX19.__captajfeoiwjri3qjr03q(list_link_cewur82214849, alpha_proteins)
            list_final_link_cewur82214849 = list_link_cewur82214849[:dp.myCFG.D11_NUMBER_TOP_RESULT]
            res.append(list_final_link_cewur82214849)
        assert (len(res) == len(spectrum_list))
        return res

    @staticmethod
    def __captainmvsijefoijrqo(search_alpha_info:CSearchFromBeta, pep_index_G, spectrum:CSpectrum, alpha_rew:list, alpha_rew_ind):
        start_time = time.time()
        alpha_mass = spectrum.mass - search_alpha_info.beta_info.mass - ATOM_MASS_P
        if search_alpha_info.dp.myCFG.C1_TYPE_TOL_PRECURSOR:
            tol_mass = spectrum.mass * search_alpha_info.dp.myCFG.C2_PPM_TOL_PRECURSOR
        else:
            tol_mass = search_alpha_info.dp.myCFG.C2_PPM_TOL_PRECURSOR
        G = alpha_mass - tol_mass
        T = alpha_mass + tol_mass
        if G > alpha_rew[-1].mass:
            return []
        if G < 0:
            G = 0.0
        if T > alpha_rew[-1].mass:
            T = alpha_rew[-1].mass
        s_start_ind = alpha_rew_ind[int((G - pep_index_G) * search_alpha_info.dp.myCFG.D12_MULTI_MASS)]
        while s_start_ind < len(alpha_rew) and alpha_rew[s_start_ind].mass < G:
            s_start_ind += 1
        s_end_ind = alpha_rew_ind[int((T - pep_index_G) * search_alpha_info.dp.myCFG.D12_MULTI_MASS)]
        while s_end_ind < len(alpha_rew) and alpha_rew[s_end_ind].mass <= T:
            s_end_ind += 1
        if s_end_ind - s_start_ind <= 0:
            return []
        candidate_peptide_list = alpha_rew[s_start_ind: s_end_ind]  # 获取候选肽段区间
        end_time = time.time()
        # print("[Debug]Time getting candidate peps for spectrum {0} is {1}".format(spectrum.title, (end_time - start_time)))
        return candidate_peptide_list

    @staticmethod
    def __capjfewij0iwji3j340apowqpokdw(alpha_cewur82214849:CPeptideInfo, alpha_proteins:CProtein, set_link_site:str, C_str:str= "", N_str:str= ""):
        alpha_peptide_sq = alpha_proteins.sq[alpha_cewur82214849.pro_index][alpha_cewur82214849.start_pos:alpha_cewur82214849.end_pos]
        alpha_vnbbmx = []
        mod_site_dict = set()
        for m in alpha_cewur82214849.mods:
            mod_site_dict.add(m.site - 1)
        for s in range(len(alpha_peptide_sq)):
            if alpha_peptide_sq[s] not in set_link_site: continue
            if s in mod_site_dict: continue
            pro_N = False
            pro_C = False
            if s == 0:
                if alpha_cewur82214849.start_pos == 0:
                    pro_N = True
            if s == len(alpha_peptide_sq) - 1:
                if alpha_cewur82214849.end_pos == alpha_cewur82214849.pro_len:
                    pro_C = True
            if not toolCheckLinkSiteValid(C_str, N_str, alpha_peptide_sq, s, pro_N, pro_C): continue
            alpha_vnbbmx.append(s)
        return alpha_vnbbmx

    @staticmethod
    def Gefewjiojfoiewjroi3wq(dp:CDataPack, spectrum:CSpectrum, alpha_cewur82214849:CPeptideInfo, alpha_link_site:int, alpha_peptide_sq:str, beta_cewur82214849:CPeptideInfo, beta_link_site:int, beta_peptide_sq:str, beta_score_feature:CSingleScoreFeature, charge_list:list=[1, 2]):
        link_cewur82214849 = CLinkPeptide(alpha_cewur82214849, alpha_link_site, beta_cewur82214849, beta_link_site)

        if beta_score_feature.uaa_clv:
            alpha_score_feature = CFunctionX2.wejwehrE34320UFJFSJOJSGOI(dp, spectrum, alpha_cewur82214849, alpha_peptide_sq, alpha_link_site, dp.myCFG.D12_MULTI_MASS, charge_list=charge_list)
            link_cewur82214849.clv = 1
        else:
            alpha_score_feature = CFunctionX2.EFHU3H2JFOIWEO3EWR32fdewfwwf(dp, spectrum, alpha_cewur82214849, alpha_peptide_sq, alpha_link_site, dp.myCFG.D12_MULTI_MASS, beta_cewur82214849.mass, charge_list=charge_list, alpha_or_beta="α")

        link_cewur82214849.score = alpha_score_feature.score + beta_score_feature.score
        link_cewur82214849.match_score = alpha_score_feature.match_score + beta_score_feature.match_score
        link_cewur82214849.alpha_score = alpha_score_feature
        link_cewur82214849.beta_score = beta_score_feature
        link_cewur82214849.alpha_sq = alpha_peptide_sq
        link_cewur82214849.beta_sq = beta_peptide_sq
        link_cewur82214849.abs_delta_da = abs(spectrum.mass - (alpha_cewur82214849.mass + beta_cewur82214849.mass + ATOM_MASS_P) )
        if alpha_cewur82214849.pro_index % 2 != 0 :
            link_cewur82214849.TD = -1
        return link_cewur82214849

    @staticmethod
    def __captajfeoiwjri3qjr03q(list_link_cewur82214849:list, alpha_proteins:CProtein):

        if len(list_link_cewur82214849) <= 1:
            return list_link_cewur82214849
        # print(pep_list[0].alpha_peptide.pro_index)
        ac = alpha_proteins.ac[list_link_cewur82214849[0].alpha_peptide.pro_index]
        if not ac.startswith("REV_"):
            return list_link_cewur82214849
        if list_link_cewur82214849[0].score - list_link_cewur82214849[1].score < 0.0001:
            return list_link_cewur82214849
        for i in range(1, len(list_link_cewur82214849)):
            if list_link_cewur82214849[i].score < list_link_cewur82214849[0].score:
                break
            new_ac = alpha_proteins.ac[0][list_link_cewur82214849[i].alpha_peptide.pro_index]
            if not new_ac.startswith("REV_"):
                break  # target protein
            if i < len(list_link_cewur82214849) and list_link_cewur82214849[i].score - list_link_cewur82214849[0].score <= 0.0001:
                tmp = list_link_cewur82214849[0]
                list_link_cewur82214849[0] = list_link_cewur82214849[i]
                list_link_cewur82214849[i] = tmp


class CFunctionX20:

    @staticmethod
    def qwe2qesadfsdfdzgwgewq(input_data_list, type_thread:int=0, number_thread:int=1):
        if len(input_data_list) < number_thread:
            number_thread = len(input_data_list)
        if type_thread == 0:
            pool = Pool(processes=number_thread)
        else:
            pool = ThreadPool(number_thread, input_data_list)
        res = pool.map(CFunctionX20.fefuoiewjiowq3jeioq, input_data_list)
        pool.close()
        pool.join()
        return res

    @staticmethod
    def fefuoiewjiowq3jeioq(search_alpha_info:CSearchFromBeta):

        spectrum_list = search_alpha_info.list_spectrum_info
        beta_cewur82214849 = search_alpha_info.beta_info
        beta_peptide_sq = search_alpha_info.beta_sq
        dp = search_alpha_info.dp
        res = []
        beta_link_site = beta_peptide_sq.find(dp.myCFG.B8_UAA_AA)
        for spectrum in spectrum_list:
            beta_score_feature = CFunctionX2.FEJWIFEehiuwfuwienfiuw38409(dp, spectrum, beta_cewur82214849, beta_peptide_sq, beta_link_site, dp.myCFG.D12_MULTI_MASS)
            open_beta_result = CFunctionX20.__captainGeJFEOIWJW(beta_cewur82214849, spectrum, beta_score_feature, beta_peptide_sq, beta_link_site)
            res.append(open_beta_result)
        return res

    @staticmethod
    def __captainGeJFEOIWJW(beta_cewur82214849:CPeptideInfo, spectrum:CSpectrum, beta_score_feature:CSingleScoreFeature, beta_peptide_sq:str, beta_link_site:int):
        return CSinglePeptide(beta_cewur82214849, beta_peptide_sq, beta_link_site, (spectrum.mass - ATOM_MASS_P - beta_cewur82214849.mass), beta_cewur82214849.pro_index, beta_score_feature.score)