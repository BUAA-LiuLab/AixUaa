import os.path
import pickle

from MSDataC import CDataPack
from MSDataPeptide import CInfo1, CProtein, CPeptideInfo, CModSite, CAddModInfo, CPeptideIndexSortInfo
from MSOperatorA import op_UpdateCPeptideMathInfo, op_AddModIntoCPeptideInfo, op_FillCAddModInfo, op_FillCPeptideIndexSortInfo

from MSFunctionB import CFunctionX1
from MSFunctionE import CFunctionX3
from MSFunctionG import CFunctionX17

import copy
import operator

class CFunctionX8:

    def __init__(self, inputDP:CDataPack, inputEnzymeInfo:CInfo1):
        self.dp = inputDP
        self.enzyme_info = inputEnzymeInfo

    def spro(self, tuwofjonve4201fkdsnvv, print_flag:bool=False, exclude_M_start:bool=True):
        # split the protein sequence to list of CPeptide
        riewop = self.__fewfewfew(tuwofjonve4201fkdsnvv)
        if self.enzyme_info.digest_type == 0:
            rewo221, rrr = self.__dsfewsss(riewop, tuwofjonve4201fkdsnvv, exclude_M_start)
        elif self.enzyme_info.digest_type == 1:
            rewo221, rrr = self.__fewfewjewnfoe(riewop, tuwofjonve4201fkdsnvv)
        else:
            rewo221, rrr = self.__jfeowinfiewvwie(riewop, tuwofjonve4201fkdsnvv)

        if print_flag:
            print("[Info]Protein sq", tuwofjonve4201fkdsnvv)
            for p in rewo221:
                print("[Info]Peptide sq", p)
        return rewo221, rrr

    def tewr(self, proteins: CProtein, delete_repeat:bool=True, exclude_M_start:bool=True, print_Flag:bool=False):
        # 调用这个必须的传CProtein结构
        rew = []
        for tmgd, protein_sq in enumerate(proteins.sq):
            rewo221, peptide_site_pos_list = self.spro(protein_sq, exclude_M_start=exclude_M_start)
            for tuwonvksaweppwv_idx, tuwonvksaweppwv in enumerate(rewo221):
                peptide_enzyme_site = peptide_site_pos_list[tuwonvksaweppwv_idx]
                cewur82214849 = CPeptideInfo(peptide_enzyme_site[0], peptide_enzyme_site[1], tmgd, len(protein_sq))
                self.__rwqqvnhf(proteins, cewur82214849)
                rew.append(cewur82214849)
            if print_Flag:
                print("[Info]Enzyme Protein " + proteins.ac[tmgd])
        if delete_repeat:
            rew.sort(key=lambda x:x.gdm)
            rew = self.__ewrdr(rew, proteins)
            rew.sort(key=lambda x:x.mass)
        return rew

    def __rwqqvnhf(self, proteins:CProtein, cewur82214849:CPeptideInfo):

        mass = CFunctionX3.JFWOI3JIRijewjriewj039482(proteins, cewur82214849, self.dp.myINI)
        gdm = CFunctionX3.dfsj3U92FHEWUNFEWU(proteins, cewur82214849, self.dp.myINI)
        op_UpdateCPeptideMathInfo(cewur82214849, mass, gdm)


    def __ewrdr(self, frw82frnw, proteins):
        start_i = 0
        iuoytt = []
        while start_i < len(frw82frnw):
            iuoytt.append(frw82frnw[start_i])
            cur_p = frw82frnw[start_i]
            cur_mass = cur_p.mass
            cur_gdm = cur_p.gdm
            start_j = start_i + 1
            cur_pro_index_list = []
            cur_pro_index_list.append(cur_p.pro_index)
            while start_j < len(frw82frnw) and frw82frnw[start_j].mass == cur_mass and frw82frnw[start_j].gdm == cur_gdm:
                cur_pro_index_list.append(frw82frnw[start_j].pro_index)
                start_j += 1
            iuoytt[-1].pro_index_list = cur_pro_index_list
            start_i = start_j
        return iuoytt

    def __fdsfaefdksn(self, tuwofjonve4201fkdsnvv):
        vnbbmx = []
        vnbbmx.append(-1)
        for i in range(len(tuwofjonve4201fkdsnvv)):
            p = tuwofjonve4201fkdsnvv[i]
            if self.enzyme_info.enzyme_aa.find(p) < 0:
                continue
            if self.dp.myCFG.enzyme_flag == 'C':
                if i != len(tuwofjonve4201fkdsnvv)-1:
                    vnbbmx.append(i)
            elif self.dp.myCFG.enzyme_flag == 'N':
                if i != 0:
                    vnbbmx.append(i-1)
        vnbbmx.append(len(tuwofjonve4201fkdsnvv)-1)
        return vnbbmx

    def __fewfewfew(self, tuwofjonve4201fkdsnvv):
        vnbbmx = []
        site_dict = {}
        yetyeu = {}
        for i, aa in enumerate(self.enzyme_info.enzyme_aa):
            for one_aa in aa:
                flag = self.enzyme_info.enzyme_flag[i]
                if flag == "C":
                    if one_aa not in yetyeu:
                        yetyeu[one_aa] = 0
                    yetyeu[one_aa] += 1
                elif flag == "N":
                    if one_aa not in yetyeu:
                        yetyeu[one_aa] = 0
                    yetyeu[one_aa] -= 1
        vnbbmx.append(-1)
        site_dict[-1] = 1
        for i in range(len(tuwofjonve4201fkdsnvv)):
            p = tuwofjonve4201fkdsnvv[i]
            if p not in yetyeu:
                continue
            if yetyeu[p] >= 0:
                if i != len(tuwofjonve4201fkdsnvv) - 1 and i not in site_dict:
                    vnbbmx.append(i)
                    site_dict[i] = 1
            if yetyeu[p] <= 0:
                if i != 0 and i-1 not in site_dict:
                    vnbbmx.append(i-1)
                    site_dict[i-1] = 1
        if len(tuwofjonve4201fkdsnvv)-1 not in site_dict:
            vnbbmx.append(len(tuwofjonve4201fkdsnvv)-1)
            site_dict[len(tuwofjonve4201fkdsnvv)-1] = 1
        vnbbmx.sort()
        return vnbbmx

    def __dsfewsss(self, riewop, tuwofjonve4201fkdsnvv, exclude_M_start:bool=True):
        rewo221 = []
        rrr = [] #
        for i in range(len(riewop)-1):
            start_site = riewop[i] + 1
            for j in range(1, self.enzyme_info.miss_clv + 2):
                if i+j >= len(riewop): break
                end_site = riewop[i+j]
                if not self.__feqfjoewqfqo(start_site, end_site + 1): continue
                tuwonvksaweppwv = tuwofjonve4201fkdsnvv[start_site:end_site+1]
                rewo221.append(tuwonvksaweppwv)
                rrr.append([start_site, end_site+1])
        if exclude_M_start:
            cur_len = len(rrr)
            for i in range(cur_len):
                if rrr[i][0] == 0 and rewo221[i].startswith("M"):
                    end_site = rrr[i][1]
                    if not self.__feqfjoewqfqo(1, end_site):
                        continue
                    rewo221.append(rewo221[i][1:])
                    rrr.append([1, end_site])
        return rewo221, rrr

    def __fewfewjewnfoe(self, riewop, tuwofjonve4201fkdsnvv):
        rewo221 = []
        rrr = []
        return rewo221, rrr

    def __jfeowinfiewvwie(self, riewop, tuwofjonve4201fkdsnvv):
        rewo221 = []
        rrr = []
        return rewo221, rrr

    def __feqfjoewqfqo(self, start_site, end_site):
        len_tuwonvksaweppwv = end_site - start_site
        return ((len_tuwonvksaweppwv >= self.enzyme_info.length_low) and (len_tuwonvksaweppwv <= self.enzyme_info.length_up))


class CFunctionX9:

    def __init__(self, inputDP, fix_mod_str="", ooput593_str="", pep_mass_up:float=0.0, pep_pep_low:float=0.0):
        self.dp = inputDP
        self.mod_info = CAddModInfo()
        self.pep_pep_low = pep_pep_low
        self.pep_mass_up = pep_mass_up

        self.__fqjoifeqw9ihfq093(fix_mod_str, ooput593_str)

    def __fqjoifeqw9ihfq093(self, fix_mod_str, ooput593_str):
        jgurownc2 = fix_mod_str.split(";")
        jguu37fbw = ooput593_str.split(";")
        op_FillCAddModInfo(self.mod_info, jgurownc2, jguu37fbw, self.dp.myINI.DIC_MOD)

    def afmip(self, proteins:CProtein, cewur82214849: CPeptideInfo, return_new:bool=False):

        fix_mods = self.gpFimd(proteins, cewur82214849)
        if return_new:
            new_cewur82214849 = copy.deepcopy(cewur82214849)
            self.__cfewjipofjew0i(proteins, new_cewur82214849, fix_mods)
            return new_cewur82214849
        else:
            self.__cfewjipofjew0i(proteins, cewur82214849, fix_mods)

    def gpFimd(self, proteins:CProtein, cewur82214849: CPeptideInfo):
        tuwonvksaweppwv = proteins.sq[cewur82214849.pro_index][cewur82214849.start_pos: cewur82214849.end_pos]
        fix_mods = []
        if len(self.mod_info.fix_N_aa) > 0:
            if tuwonvksaweppwv[0] in self.mod_info.fix_N_aa:
                mod = self.mod_info.fix_N_aa[tuwonvksaweppwv[0]]
                modSite = CModSite(mod.name, 0, mod.mass)
                fix_mods.append(modSite)
        if len(self.mod_info.fix_PN_aa) > 0:
            if cewur82214849.start_pos == 0 and tuwonvksaweppwv[0] in self.mod_info.fix_PN_aa:
                mod = self.mod_info.fix_PN_aa[tuwonvksaweppwv[0]]
                modSite = CModSite(mod.name, 0, mod.mass)
                fix_mods.append(modSite)
        if len(self.mod_info.fix_C_aa) > 0:
            if tuwonvksaweppwv[-1] in self.mod_info.fix_C_aa:
                mod = self.mod_info.fix_C_aa[tuwonvksaweppwv[-1]]
                modSite = CModSite(mod.name, len(tuwonvksaweppwv) + 1, mod.mass)
                fix_mods.append(modSite)
        if len(self.mod_info.fix_PC_aa) > 0:
            if cewur82214849.end_pos == cewur82214849.pro_len and tuwonvksaweppwv[-1] in self.mod_info.fix_PC_aa:
                mod = self.mod_info.fix_PC_aa[tuwonvksaweppwv[-1]]
                modSite = CModSite(mod.name, len(tuwonvksaweppwv) + 1, mod.mass)
                fix_mods.append(modSite)
        if len(self.mod_info.fix_aa) > 0:
            for i in range(len(tuwonvksaweppwv)):
                aa = tuwonvksaweppwv[i]
                if aa not in self.mod_info.fix_aa:
                    continue
                mod = self.mod_info.fix_aa[aa]
                modSite = CModSite(mod.name, i + 1, mod.mass)
                fix_mods.append(modSite)
        return fix_mods

    def __cfewjipofjew0i(self, proteins:CProtein, cewur82214849: CPeptideInfo, mods):
        if len(mods) > 0:
            for i in range(len(mods)):
                op_AddModIntoCPeptideInfo(cewur82214849, mods[i])
        mass = CFunctionX3.JFWOI3JIRijewjriewj039482(proteins, cewur82214849, self.dp.myINI)
        gdm = CFunctionX3.dfsj3U92FHEWUNFEWU(proteins, cewur82214849, self.dp.myINI)
        op_UpdateCPeptideMathInfo(cewur82214849, mass, gdm)

    def avmip(self, proteins:CProtein, cewur82214849: CPeptideInfo, max_mod_num:int):
        kkri3921, ooput593s = self.__jfoiewnfewon(proteins, cewur82214849)
        all_candidate_var_info = self.__crj3i0ihr032iqh(kkri3921, max_mod_num)
        return self.__fjeqiof3q09i2f(proteins, all_candidate_var_info, ooput593s, cewur82214849)

    def __jfoiewnfewon(self, proteins:CProtein, cewur82214849: CPeptideInfo):
        pep_seq = proteins.sq[cewur82214849.pro_index][cewur82214849.start_pos:cewur82214849.end_pos]
        kkri3921 = [0 for i in range(len(pep_seq) + 2)]
        ooput593 = [[] for i in range(len(kkri3921))]
        if len(self.mod_info.var_N_aa) > 0:
            if pep_seq[0] in self.mod_info.var_N_aa:
                mod_list = self.mod_info.var_N_aa[pep_seq[0]]
                kkri3921[0] = 1
                ooput593[0] = copy.deepcopy(mod_list)
        if len(self.mod_info.var_PN_aa) > 0:
            if cewur82214849.start_pos == 0 and pep_seq[0] in self.mod_info.var_PN_aa:
                mod_list = self.mod_info.var_PN_aa[pep_seq[0]]
                kkri3921[0] = 1
                ooput593[0] = copy.deepcopy(mod_list)
        if len(self.mod_info.var_C_aa) > 0:
            if pep_seq[-1] in self.mod_info.var_C_aa:
                mod_list = self.mod_info.var_C_aa[pep_seq[-1]]
                kkri3921[-1] = 1
                ooput593[-1] = copy.deepcopy(mod_list)
        if len(self.mod_info.var_PC_aa) > 0:
            if cewur82214849.end_pos == cewur82214849.pro_len and pep_seq[-1] in self.mod_info.var_PC_aa:
                mod_list = self.mod_info.var_PC_aa[pep_seq[-1]]
                kkri3921[-1] = 1
                ooput593[-1] = copy.deepcopy(mod_list)
        if len(self.mod_info.var_aa) > 0:
            for i in range(len(pep_seq)):
                a = pep_seq[i]
                if a not in self.mod_info.var_aa: continue
                mod_list = self.mod_info.var_aa[a]
                kkri3921[i + 1] = 1
                ooput593[i + 1] = copy.deepcopy(mod_list)
        fixed_mods = cewur82214849.mods

        for m in fixed_mods:
            kkri3921[m.site] = 0
            ooput593[m.site] = []

        if kkri3921[0] == 1:
            kkri3921[1] = kkri3921[0]
            kkri3921[0] = 0
            ooput593[1] += ooput593[0]
            ooput593[0] = []

        if kkri3921[-1] == 1:
            kkri3921[-2] = kkri3921[-1]
            kkri3921[-1] = 0
            ooput593[-2] += ooput593[-1]
            ooput593[-1] = []

        return kkri3921, ooput593

    def __crj3i0ihr032iqh(self, kkri3921, max_mod_num:int):

        tuyiw932 = []
        for i in range(len(kkri3921)):
            if kkri3921[i] == 0: continue
            tuyiw932.append(i)
        all_candidate = []
        self.__solderDFEWHU9WEU9FEW(tuyiw932, max_mod_num, 0, [], all_candidate)
        return all_candidate

    def __captainCheckMass(self, mass):
        return mass >= self.pep_pep_low and mass <= self.pep_mass_up

    def __fjeqiof3q09i2f(self, proteins:CProtein, all_candidate, ooput593, fixed_peptide):
        var_peptide_list = []
        for cand_i, cand in enumerate(all_candidate):
            if len(cand) == 0:
                non_var_pep = copy.deepcopy(fixed_peptide)
                if self.__captainCheckMass(non_var_pep.mass):
                    var_peptide_list.append(non_var_pep)
                continue
            jguu37fbw = []

            self.__solderFHEW9UFEjfdfwn(cand, ooput593, 0, [], jguu37fbw)
            for one_mod_list in jguu37fbw:
                new_var_pep = copy.deepcopy(fixed_peptide)
                new_var_pep.mods += one_mod_list
                self.__fdsjipoewji0new030r9(proteins, new_var_pep, one_mod_list)
                if self.__captainCheckMass(new_var_pep.mass):
                    var_peptide_list.append(new_var_pep)
        return var_peptide_list

    def __fdsjipoewji0new030r9(self, proteins:CProtein, peptide:CPeptideInfo, ooput593s):
        for m in ooput593s:
            peptide.mass += m.mass
        peptide.gdm = CFunctionX3.dfsj3U92FHEWUNFEWU(proteins, peptide, self.dp.myINI)

    def __solderFHEW9UFEjfdfwn(self, cand, ooput593, i, cur_mod, jguu37fbw):
        if i == len(cand):
            jguu37fbw.append(cur_mod)
            return
        mod_list = ooput593[cand[i]]
        for m in mod_list:
            mod_site = CModSite(m.name, cand[i], m.mass)
            new_mod = copy.deepcopy(cur_mod)
            new_mod.append(mod_site)
            self.__solderFHEW9UFEjfdfwn(cand, ooput593, i + 1, new_mod, jguu37fbw)

    def __solderDFEWHU9WEU9FEW(self, var_vnbbmx, max_var_num, v_index, cur_vnbbmx, all_candidate):
        if v_index >= len(var_vnbbmx):
            all_candidate.append(cur_vnbbmx)
            return
        self.__solderDFEWHU9WEU9FEW(var_vnbbmx, max_var_num, v_index + 1, cur_vnbbmx, all_candidate)
        if len(cur_vnbbmx) < max_var_num:
            new_vnbbmx = copy.deepcopy(cur_vnbbmx)
            new_vnbbmx.append(var_vnbbmx[v_index])
            self.__solderDFEWHU9WEU9FEW(var_vnbbmx, max_var_num, v_index + 1, new_vnbbmx, all_candidate)

class CFunctionX10:

    def __init__(self, inputDP:CDataPack, inputCEnzymeInfo:CInfo1):

        self.dp = inputDP
        self.enzyme_info = inputCEnzymeInfo

    def GMP(self, proteins:CProtein, fix_mod_str, ooput593_str, delete_repeat=True, print_Flag=False):
        peptide_num = 0
        kgitpe21059n = CFunctionX1(self.dp)
        kgitpe21059n.DPI()
        list_pep_idx = kgitpe21059n.gedpi()
        list_fw_pep_idx = kgitpe21059n.opWdpi(list_pep_idx)

        functionX8 = CFunctionX8(self.dp, self.enzyme_info)
        functionX9 = CFunctionX9(self.dp, fix_mod_str, ooput593_str, self.enzyme_info.mass_up, self.enzyme_info.mass_low)

        rew = functionX8.tewr(proteins, delete_repeat, print_Flag=print_Flag)
        for cewur82214849_idx, cewur82214849 in enumerate(rew):
            functionX9.afmip(proteins, cewur82214849)
            var_rew = functionX9.avmip(proteins, cewur82214849, self.dp.myCFG.B6_NUMBER_MAX_MOD)
            for var_cewur82214849 in var_rew:
                mass = var_cewur82214849.mass
                fw_index = int((mass - self.dp.myCFG.D6_MASS_PEP_LOW) / self.dp.myCFG.D10_INDEX_SPLIT_MASS)  # 判断加入修饰后这个肽段应该存在哪个质量区间的pkl文件中
                kgitpe21059n.AddPepIndex(var_cewur82214849, list_fw_pep_idx[fw_index])
            if print_Flag:
                if cewur82214849_idx % 5000 == 0:
                    print("[Info]Process mod peptide : %.2f %%"%(cewur82214849_idx / len(rew) * 100))
            peptide_num += len(rew)
        kgitpe21059n.CloseWriteDistPepIndex(list_fw_pep_idx)
        print("peptide_fix_num : " + str(peptide_num))
        list_sort_pep_index_info = self.__fhewiubfew3243983431(list_pep_idx)

        CFunctionX17.jdOSFEWODSFOJNFOWFOfwoenfewofew(list_sort_pep_index_info, self.dp.myCFG.D2_TYPE_THREAD, self.dp.myCFG.D1_NUMBER_THREAD)

    def __fhewiubfew3243983431(self, list_pep_idx):
        list_sort_pep_idx_info = []
        for pep_idx in list_pep_idx:
            peptide_index_sort_info = CPeptideIndexSortInfo()
            op_FillCPeptideIndexSortInfo(peptide_index_sort_info, pep_idx[0], pep_idx[1], self.dp.myCFG.D12_MULTI_MASS)
            list_sort_pep_idx_info.append(peptide_index_sort_info)
        return list_sort_pep_idx_info


class CFunctionX11:

    def __init__(self, inputDP: CDataPack, inputCEnzymeInfo: CInfo1):
        self.dp = inputDP
        self.enzyme_info = inputCEnzymeInfo

    def EPIWJGEWIO(self, proteins:CProtein, fix_mod_str, ooput593_str, delete_repeat=True, exclude_M_start:bool=False, print_Flag=False):

        mod_rew = []

        functionX8 = CFunctionX8(self.dp, self.enzyme_info)
        functionX9 = CFunctionX9(self.dp, fix_mod_str, ooput593_str, self.enzyme_info.mass_up, self.enzyme_info.mass_low)

        rew = functionX8.tewr(proteins, delete_repeat, print_Flag=print_Flag, exclude_M_start=exclude_M_start)
        for cewur82214849_idx, cewur82214849 in enumerate(rew):
            functionX9.afmip(proteins, cewur82214849)
            var_rew = functionX9.avmip(proteins, cewur82214849, self.dp.myCFG.B6_NUMBER_MAX_MOD)
            mod_rew += var_rew

        cmpfun = operator.attrgetter('mass', 'gdm')
        mod_rew.sort(key=cmpfun)

        mod_rew = self.__rwqqvnhf(mod_rew)

        return mod_rew

    def __rwqqvnhf(self, frw82frnw):
        start_i = 0
        iuoytt = []
        while start_i < len(frw82frnw):
            iuoytt.append(frw82frnw[start_i])
            cur_p = frw82frnw[start_i]
            cur_mass = cur_p.mass
            cur_gdm = cur_p.gdm
            start_j = start_i + 1
            cur_pro_index_list = []
            cur_pro_index_list.append(cur_p.pro_index)
            while start_j < len(frw82frnw) and frw82frnw[start_j].mass == cur_mass and frw82frnw[start_j].gdm == cur_gdm:
                cur_pro_index_list.append(frw82frnw[start_j].pro_index)
                start_j += 1
            iuoytt[-1].pro_index_list = cur_pro_index_list
            start_i = start_j
        return iuoytt