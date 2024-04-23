from MSLogging import CLogging
from MSSystem import MOLECULE_MASS_H2O, AA_NUM, ATOM_MASS_P, AVRG_AA_MASS
from MSDataC import CDataPack, CINI
from MSDataPeptide import CModSite, CMod, CInfo1, CProtein, CPeptideInfo, CAddModInfo
from MSDataD import CSpectrum
from MSDataB import CSingleScoreFeature

import numpy as np
import copy
import re
import math

class CFunctionX2:

    @staticmethod
    def EFHU3H2JFOIWEO3EWR32fdewfwwf(inputDP:CDataPack, spectrum:CSpectrum, cewur82214849:CPeptideInfo, peptide_sq:str, link_site:int, multi:int=1, link_mass=None, charge_list:list=[1, 2], alpha_or_beta:str= "β")->CSingleScoreFeature:
        if link_mass is None:  # mass1 is the mass of alpha peptide
            link_mass = spectrum.mass - ATOM_MASS_P - cewur82214849.mass

        jjtu3922, all_peak_inten = CFunctionX2.__cajfidsnfeiw3984U2JF9EHWF9E(spectrum, multi)  # create the peak index  # 一张谱图生成一个峰强索引和计算强度总和
        beta_ion_list, beta_ion_flag_list = CFunctionX2.FWHEUFWU389y439hfuwh8u3he3e94y7(peptide_sq, cewur82214849, link_site, link_mass, inputDP.myINI.DIC_AA, charge_list, alpha_or_beta=alpha_or_beta)
        sq_len = float(len(peptide_sq))
        avrgd = cewur82214849.mass / AVRG_AA_MASS
        norm = (sq_len * sq_len) / (avrgd * avrgd)

        single_score = CFunctionX2.__MI4i3u4320jfiewniri324j3(jjtu3922, spectrum.peaks, all_peak_inten,
                                                              beta_ion_list, beta_ion_flag_list, inputDP.myCFG.C3_TYPE_TOL_FRAGMENT,
                                                              inputDP.myCFG.C4_PPM_TOL_FRAGMENT, multi, norm, charge_list=charge_list)
        return single_score

    @staticmethod
    def EFHIUWHFU32Nfuewinwe2134(inputDP:CDataPack, spectrum:CSpectrum, cewur82214849:CPeptideInfo, peptide_sq:str, link_site:int, multi:int=1)->CSingleScoreFeature:

        jjtu3922, all_peak_inten = CFunctionX2.__cajfidsnfeiw3984U2JF9EHWF9E(spectrum, multi)  # create the peak index  # 一张谱图生成一个峰强索引和计算强度总和

        llkg0922, beta_clv_ion_flag_list = CFunctionX2.iwuhuhfu9h9u32h4u32h4f(peptide_sq, cewur82214849, link_site, inputDP.myINI.CLV_UAA_MASS, inputDP.myINI.DIC_AA)

        sq_len = float(len(peptide_sq))
        avrgd = cewur82214849.mass / AVRG_AA_MASS
        norm = (sq_len * sq_len) / (avrgd * avrgd)

        single_clv_score = CFunctionX2.__MI4i3u4320jfiewniri324j3(jjtu3922, spectrum.peaks, all_peak_inten, llkg0922, beta_clv_ion_flag_list, inputDP.myCFG.C3_TYPE_TOL_FRAGMENT, inputDP.myCFG.C4_PPM_TOL_FRAGMENT, multi, norm, uaa_clv=True)
        return single_clv_score

    @staticmethod
    def wejwehrE34320UFJFSJOJSGOI(inputDP:CDataPack, spectrum:CSpectrum, cewur82214849:CPeptideInfo, peptide_sq:str, link_site:int, multi:int=1, charge_list:list=[1, 2])->CSingleScoreFeature:

        jjtu3922, all_peak_inten = CFunctionX2.__cajfidsnfeiw3984U2JF9EHWF9E(spectrum, multi)  # create the peak index  # 一张谱图生成一个峰强索引和计算强度总和

        powrj22, orwpjwvnfnbei32 = CFunctionX2.FWHEUFWU389y439hfuwh8u3he3e94y7(peptide_sq, cewur82214849, link_site, 0.0, inputDP.myINI.DIC_AA, alpha_or_beta="α", charge_list=charge_list)

        sq_len = float(len(peptide_sq))
        avrgd = cewur82214849.mass / AVRG_AA_MASS
        norm = (sq_len * sq_len) / (avrgd * avrgd)

        single_clv_score = CFunctionX2.__MI4i3u4320jfiewniri324j3(jjtu3922, spectrum.peaks, all_peak_inten, powrj22, orwpjwvnfnbei32, inputDP.myCFG.C3_TYPE_TOL_FRAGMENT, inputDP.myCFG.C4_PPM_TOL_FRAGMENT, multi, norm, uaa_clv=True)
        return single_clv_score

    @staticmethod
    def FEJWIFEehiuwfuwienfiuw38409(inputDP: CDataPack, spectrum: CSpectrum, cewur82214849: CPeptideInfo, peptide_sq: str, link_site: int, multi: int = 1, link_mass=None) -> CSingleScoreFeature:
        if link_mass is None:  # mass1 is the mass of alpha peptide
            link_mass = spectrum.mass - ATOM_MASS_P - cewur82214849.mass

        jjtu3922, all_peak_inten = CFunctionX2.__cajfidsnfeiw3984U2JF9EHWF9E(spectrum,
                                                                             multi)  # create the peak index  # 一张谱图生成一个峰强索引和计算强度总和
        beta_ion_list, beta_ion_flag_list = CFunctionX2.FWHEUFWU389y439hfuwh8u3he3e94y7(peptide_sq, cewur82214849, link_site, link_mass, inputDP.myINI.DIC_AA)

        sq_len = float(len(peptide_sq))
        avrgd = cewur82214849.mass / AVRG_AA_MASS
        norm = (sq_len * sq_len) / (avrgd * avrgd)

        single_score = CFunctionX2.__MI4i3u4320jfiewniri324j3(jjtu3922, spectrum.peaks,
                                                              all_peak_inten, beta_ion_list,
                                                              beta_ion_flag_list,
                                                              inputDP.myCFG.C3_TYPE_TOL_FRAGMENT,
                                                              inputDP.myCFG.C4_PPM_TOL_FRAGMENT, multi,
                                                              norm)
        return single_score

    @staticmethod
    def FWHEUFWU389y439hfuwh8u3he3e94y7(pep_sq: str, pep_info: CPeptideInfo, link_site: int, add_mass: float, aa2mass: dict, charge_list:list=[1, 2], alpha_or_beta:str= "β"):
        aa2 = CFunctionX2.vonjdsoiruiui32u43(pep_sq, pep_info.mods, aa2mass)

        b2_mz, y2_mz = [], []
        b2_flag, y2_flag = [], []
        tmp_mass = 0.0
        for i in range(len(pep_sq) - 1):
            tmp_mass += aa2[i]
            if i < link_site:
                b_charge0 = tmp_mass
            else:
                b_charge0 = tmp_mass + add_mass
            for charge in charge_list:
                b_charge = b_charge0 + ATOM_MASS_P * charge
                b2_mz.append(b_charge / charge)
                b2_flag.append(alpha_or_beta + "b" + str(i + 1) + "+" * charge)
        tmp_mass = 0.0
        for i in range(1, len(pep_sq))[::-1]:
            tmp_mass += aa2[i]
            if i > link_site:
                y_charge0 = tmp_mass + MOLECULE_MASS_H2O
            else:
                y_charge0 = tmp_mass + add_mass + MOLECULE_MASS_H2O
            for charge in charge_list:
                y_charge = y_charge0 + (ATOM_MASS_P * charge)
                y2_mz.append(y_charge / charge)
                y_pos = len(pep_sq) - i
                y2_flag.append(alpha_or_beta + "y" + str(y_pos) + "+" * charge)
        return b2_mz + y2_mz, b2_flag + y2_flag

    @staticmethod
    def iwuhuhfu9h9u32h4u32h4f(pep_sq:str, pep_info:CPeptideInfo, link_site: int, uaa_mass: float, aa2mass:dict):
        aa2 = CFunctionX2.vonjdsoiruiui32u43(pep_sq, pep_info.mods, aa2mass)

        b2_mz, y2_mz = [], []
        b2_flag, y2_flag = [], []
        tmp_mass = 0.0
        for i in range(len(pep_sq) - 1):
            if i == link_site:
                tmp_mass += uaa_mass
            else:
                tmp_mass += aa2[i]
            b_charge1 = tmp_mass + ATOM_MASS_P
            b_charge2 = (b_charge1 + ATOM_MASS_P) / 2
            b2_mz.append(b_charge1)
            b2_mz.append(b_charge2)
            b2_flag.append("βb" + str(i + 1) + "+")
            b2_flag.append("βb" + str(i + 1) + "++")
        tmp_mass = 0.0
        for i in range(1, len(pep_sq))[::-1]:
            if i == link_site:
                tmp_mass += uaa_mass
            else:
                tmp_mass += aa2[i]
            y_charge1 = tmp_mass + ATOM_MASS_P + MOLECULE_MASS_H2O
            y_charge2 = (y_charge1 + ATOM_MASS_P) / 2
            y2_mz.append(y_charge1)
            y2_mz.append(y_charge2)
            y_pos = len(pep_sq) - i
            y2_flag.append("βy" + str(y_pos) + "+")
            y2_flag.append("βy" + str(y_pos) + "++")
        return b2_mz + y2_mz, b2_flag + y2_flag

    @staticmethod
    def vonjdsoiruiui32u43(sq:str, mods:list, aa2mass:dict):
        aa = [0.0 for i in range(len(sq))]
        for m in mods:
            aa[m.site-1] += m.mass
        for i, p in enumerate(sq):
            if p <'A' or p > 'Z':
                continue
            aa[i] += aa2mass[int(ord(p)-ord('A'))]
        return aa

    @staticmethod
    def __MI4i3u4320jfiewniri324j3(peak_index, peaks, all_peak_inten, ions, ion_flags, is_ppm_fra, tol_fra, multi, norm=1.0, need_compute_conscore:bool=True, uaa_clv:bool=False, charge_list:list = [1, 2], ):
        match_num = 0
        match_inten = 0.0
        match_peak_index, match_peak_me = [], []
        ion_num = len(charge_list) * 2
        match_N_flag = [0 for _ in range(int(len(ions) / ion_num))]
        match_C_flag = [0 for _ in range(int(len(ions) / ion_num))]
        match_score = 0.0
        max_me = 20.0
        if is_ppm_fra:
            max_me = tol_fra * 1e6
        else:
            max_me = tol_fra
        for mi, mz in enumerate(ions):
            if is_ppm_fra:
                tol_mass = mz * tol_fra
            else:
                tol_mass = tol_fra
            G = mz - tol_mass
            T = mz + tol_mass
            if T <= 0: continue
            if G > peaks[-1].mz: continue
            if G < 0: G = 0.0
            if T > peaks[-1].mz: T = peaks[-1].mz
            s_start_ind = peak_index[int(G * multi)]
            while s_start_ind < len(peaks) and peaks[s_start_ind].mz < G:
                s_start_ind += 1
            if int(T * multi) >= len(peak_index):
                s_end_ind = peak_index[-1]
            else:
                s_end_ind = peak_index[int(T * multi)]
            while s_end_ind < len(peaks) and peaks[s_end_ind].mz <= T:
                s_end_ind += 1
            if s_start_ind >= s_end_ind: continue
            match_num += 1
            max_inten = 0.0
            max_inten_index = -1
            for i in range(s_start_ind, s_end_ind):
                if peaks[i].inten > max_inten:
                    max_inten = peaks[i].inten
                    max_inten_index = i
            # compute Score

            if is_ppm_fra:
                one_me = (peaks[max_inten_index].mz - mz) * 1e6 / mz
            else:
                one_me = peaks[max_inten_index].mz - mz

            abs_one_me = abs(one_me)
            one_score = CFunctionX2.__JUEHR8348fsdhfuihweiuhr(norm, math.log(max_inten)) * math.cos(abs_one_me / max_me * 3.1415926 / 2)
            match_score += one_score
            match_inten += max_inten
            match_peak_index.append((max_inten_index, ion_flags[mi]))
            match_peak_me.append(one_me)
            if mi < int(len(ions) / 2):
                if mi % len(charge_list) == 0: match_N_flag[int(mi / len(charge_list))] = one_score
            else:
                if mi % len(charge_list) == 0: match_C_flag[int((len(ions) - mi) / len(charge_list)) - 1] = one_score
        match_ion_num_ratio = match_num * 1.0 / len(ions)
        match_inten_ratio = match_inten / all_peak_inten
        mean_frag_me, std_frag_me = max_me, 100.0
        if len(match_peak_me) > 0: mean_frag_me = abs(np.mean(match_peak_me))
        if len(match_peak_me) > 1: std_frag_me = np.std(match_peak_me)
        b_s, y_s, by_s = 0.0, 0.0, 0.0
        if need_compute_conscore:
            b_s, y_s, by_s = CFunctionX2.__FJW9E3u432894hfuehwn(match_N_flag, match_C_flag)
        single_match_score_feature = CSingleScoreFeature(match_score, mean_frag_me, std_frag_me, match_ion_num_ratio,
                                             match_inten_ratio,
                                             b_s, y_s, by_s, match_peak_index, uaa_clv)
        return single_match_score_feature

    @staticmethod
    def __FJW9E3u432894hfuehwn(match_N_flag, match_C_flag):
        b_s, y_s, by_s = 0, 0, 0
        TAG_LEN = 4
        for i in range(len(match_N_flag) - TAG_LEN + 1):
            bCont = True
            for j in range(TAG_LEN):
                if match_N_flag[i + j] == 0.0:
                    bCont = False
                    break
            if bCont:
                lfSum = 0.0
                for j in range(TAG_LEN): lfSum += match_N_flag[i + j]
                b_s += (lfSum / TAG_LEN)

            bCont = True
            for j in range(TAG_LEN):
                if match_C_flag[i + j] == 0.0:
                    bCont = False
                    break
            if bCont:
                lfSum = 0.0
                for j in range(TAG_LEN): lfSum += match_C_flag[i + j]
                y_s += (lfSum / TAG_LEN)

            bCont = True
            for j in range(TAG_LEN):
                if match_N_flag[i + j] == 0.0 and match_C_flag[i + j] == 0.0:
                    bCont = False
                    break
            if bCont:
                lfSum = 0.0
                for j in range(TAG_LEN): lfSum += (match_N_flag[i + j] + match_C_flag[i + j])
                by_s += (lfSum / TAG_LEN)
        return b_s, y_s, by_s

    @staticmethod
    def __JUEHR8348fsdhfuihweiuhr(len_norm, f, K1=0.08, constant_b=0.95):
        return (K1 + 1) * f / (K1 * (1 - constant_b + constant_b * len_norm) + f)

    @staticmethod
    def __cajfidsnfeiw3984U2JF9EHWF9E(spectrum:CSpectrum, multi:int):
        T = spectrum.peaks[-1].mz
        num = int(T * multi) + 1
        jjtu3922 = [-1 for i in range(num)]
        all_peak_inten = 0.0
        for i, peak in enumerate(spectrum.peaks):
            all_peak_inten += peak.inten
            index = int(peak.mz * multi)
            if index >= num: index = num - 1
            if jjtu3922[index] == -1:
                jjtu3922[index] = i
        end_val = len(spectrum.peaks)
        for i in range(num)[::-1]:
            if jjtu3922[i] == -1:
                jjtu3922[i] = end_val
            else:
                end_val = jjtu3922[i]
        return jjtu3922, all_peak_inten

class CFunctionX3:
    @staticmethod
    def soidjeiwjr0i32r3(elem_str, ini_info:CINI):
        tmp = re.split('\(|\)', elem_str)
        mass = 0.0
        for i in range(len(tmp)):
            if i % 2 != 0:
                continue
            if tmp[i].strip() == "":
                continue
            if tmp[i] not in ini_info.DIC_ELEMENT_MASS.keys():
                continue
            mass += (float(ini_info.DIC_ELEMENT_MASS[tmp[i]]) * int(tmp[i + 1]))
        return mass

    @staticmethod
    def JFWOI3JIRijewjriewj039482(proteins:CProtein, peptide:CPeptideInfo, ini_info:CINI, add_mass:float=0.0):
        tuwonvksaweppwv = proteins.sq[peptide.pro_index][peptide.start_pos: peptide.end_pos]
        mass = add_mass + MOLECULE_MASS_H2O
        for aa in tuwonvksaweppwv:
            if aa < 'A' or aa > 'Z':
                continue
            else:
                mass += ini_info.DIC_AA[int(ord(aa) - ord('A'))]
        for m in peptide.mods:
            mass += m.mass
        return mass

    @staticmethod
    def dfsj3U92FHEWUNFEWU(proteins:CProtein, cewur82214849: CPeptideInfo, ini_info:CINI):
        tuwonvksaweppwv = proteins.sq[cewur82214849.pro_index][cewur82214849.start_pos:cewur82214849.end_pos]
        aa_num = AA_NUM
        ooput593s_dict = {}
        for m in cewur82214849.mods:
            if m.mod_name not in ini_info.DIC_SET_MOD:
                CLogging.LogError(m.mod_name)
            ooput593s_dict[m.site] = ini_info.DIC_SET_MOD[m.mod_name]
        new_sq_list = []
        for i in range(len(tuwonvksaweppwv)):
            # 遍历这个肽段的氨基酸
            if i in ooput593s_dict:
                # 这个位点发生修饰
                new_sq_list.append(aa_num + ooput593s_dict[i])
            new_sq_list.append(int(ord(tuwonvksaweppwv[i]) - ord('A')))
        gdm_value = 0.0
        for i in range(len(new_sq_list)):
            ind = i
            if ind >= len(ini_info.MATRIX_GDM[0]): ind = len(ini_info.MATRIX_GDM[0]) - 1
            ind2 = new_sq_list[i]
            if ind2 >= len(ini_info.MATRIX_GDM): ind2 = len(ini_info.MATRIX_GDM) - 1
            gdm_value += ini_info.MATRIX_GDM[ind2][ind]
        return gdm_value

class CFunctionX4:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def do(self):
        self.__captainLoadElementINI()
        self.__captainLoadAAFile()
        self.__captainLoadModFile()
        self.__captainLoadModSet()
        self.__cafewh9uh329u9rjf39()
        self.__captainUAACLVMass()
        self.__captainUAALinkerMass()

    def __captainLoadElementINI(self):
        with open(self.dp.myCFG.F1_PATH_INI_ELEMENT, 'r') as f:
            for line in f:
                if line.find('=') < 0:
                    continue
                if line.strip() == "" or line[0] == '@':
                    continue
                line = line.strip()
                line = line[line.find('=')+1:]
                tmp = line.split('|')
                element_name = tmp[0]
                element_mass = [float(v) for v in tmp[1][:-1].split(',')]
                element_percentage = [float(v) for v in tmp[2][:-1].split(',')]
                assert(len(element_mass) == len(element_percentage))
                index = element_percentage.index(max(element_percentage))
                self.dp.myINI.DIC_ELEMENT_MASS[element_name] = element_mass[index]

    def __captainLoadAAFile(self):
        with open(self.dp.myCFG.F2_PATH_INI_AA, 'r') as f:
            for line in f:
                if line.find('=') < 0:
                    continue
                if line.strip() == "" or line[0] == '@':
                    continue
                name, value = line.strip().split('=')
                aa, elem_str = value[:-1].split('|')
                if aa[0] == self.dp.myCFG.B8_UAA_AA:
                    elem_str = self.dp.myCFG.B15_UAA_COM
                self.dp.myINI.DIC_AA[int(ord(aa[0]) - ord('A'))] = CFunctionX3.soidjeiwjr0i32r3(elem_str, self.dp.myINI)

    def __captainLoadModFile(self):
        with open(self.dp.myCFG.F3_PATH_INI_MOD, 'r') as f:
            for line in f:
                if line.find('=') < 0: continue
                if line.strip() == "" or line[0] == '@': continue
                if line.startswith("name"): continue
                name, value = line.strip().split('=')
                tmp = value.split(' ')
                mass = float(tmp[2])
                aa = tmp[0]
                type = tmp[1]
                self.dp.myINI.DIC_MOD[name] = CMod(name, mass, type, aa)

    def __captainLoadModSet(self):
        for alpha_fix_mod in self.dp.myCFG.B4_NAME_MOD_FIX.split(";"):
            if len(alpha_fix_mod) == 0:
                continue
            if alpha_fix_mod not in self.dp.myINI.DIC_MOD.keys():
                CLogging.LogError(alpha_fix_mod)
            if alpha_fix_mod in self.dp.myINI.DIC_SET_MOD.keys() and len(alpha_fix_mod) > 0:
                continue
            self.dp.myINI.DIC_SET_MOD[alpha_fix_mod] = len(self.dp.myINI.DIC_SET_MOD.keys())
        for alpha_ooput593 in self.dp.myCFG.B5_NAME_MOD_VAR.split(";"):
            if len(alpha_ooput593) == 0:
                continue
            if alpha_ooput593 not in self.dp.myINI.DIC_MOD.keys() and len(alpha_ooput593) > 0:
                CLogging.LogError(alpha_ooput593)
            if alpha_ooput593 in self.dp.myINI.DIC_SET_MOD.keys():
                continue
            self.dp.myINI.DIC_SET_MOD[alpha_ooput593] = len(self.dp.myINI.DIC_SET_MOD.keys())

        for beta_fix_mod in self.dp.myCFG.B13_UAA_NAME_MOD_FIX.split(";"):
            if len(beta_fix_mod) == 0:
                continue
            if beta_fix_mod not in self.dp.myINI.DIC_MOD.keys() and len(beta_fix_mod) > 0:
                CLogging.LogError(beta_fix_mod)
            if beta_fix_mod in self.dp.myINI.DIC_SET_MOD.keys():
                continue
            self.dp.myINI.DIC_SET_MOD[beta_fix_mod] = len(self.dp.myINI.DIC_SET_MOD.keys())
        for beta_ooput593 in self.dp.myCFG.B13_UAA_NAME_MOD_FIX.split(";"):
            if len(beta_ooput593) == 0:
                continue
            if beta_ooput593 not in self.dp.myINI.DIC_MOD.keys() and len(beta_ooput593) > 0:
                CLogging.LogError(beta_ooput593)
            if beta_ooput593 in self.dp.myINI.DIC_SET_MOD.keys():
                continue
            self.dp.myINI.DIC_SET_MOD[beta_ooput593] = len(self.dp.myINI.DIC_SET_MOD.keys())

    def __cafewh9uh329u9rjf39(self):
        ooput593 = self.dp.myINI.DIC_SET_MOD
        prime_num = self.dp.myCFG.D9_LEN_PEP_UP + self.dp.myCFG.B6_NUMBER_MAX_MOD + 1
        aa_num = AA_NUM + len(ooput593) + 1
        for i in range(aa_num):
            tmp = [0.0 for j in range(prime_num)]
            self.dp.myINI.MATRIX_GDM.append(tmp)
        prime_list = [0 for j in range(prime_num)]
        prime_list[0] = 2
        for i in range(1, prime_num):
            start_v = prime_list[i - 1] + 1
            while True:
                is_prime = True
                v = int(math.sqrt(start_v))
                if v == 1: v = 2
                for k in range(v, start_v):
                    if start_v % k == 0:
                        is_prime = False
                        break
                if is_prime: break
                start_v += 1
            prime_list[i] = start_v

        for i in range(aa_num):
            for j in range(prime_num):
                self.dp.myINI.MATRIX_GDM[i][j] = (i + 1) * 1.0 * math.log(prime_list[j])

    def __captainUAACLVMass(self):

        self.dp.myINI.CLV_UAA_MASS = CFunctionX3.soidjeiwjr0i32r3(self.dp.myCFG.G2_BETA_UAA_NEW_COM, self.dp.myINI)

    def __captainUAALinkerMass(self):
        self.dp.myINI.UAA_LINKERS_MASS = CFunctionX3.soidjeiwjr0i32r3(self.dp.myCFG.B15_UAA_COM, self.dp.myINI)