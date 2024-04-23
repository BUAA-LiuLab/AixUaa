from MSDataC import CConfig, CDataPack
from MSDataPeptide import CProtein
from MSDataD import CPeak, CSpectrum
from MSDataB import CLinkPeptide

from MSSystem import ATOM_MASS_P,NAME_PROTEIN_PKL
from MSLogging import CLogging

import os

try:
    import cPickle as pickle
except:
    import pickle
import operator
import copy

class CFunctionX12:

    @staticmethod
    def LoadP(path, using_decoy = True):
        pl = CProtein()
        ac, de, sq = "", "", ""
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line[0] == '>':
                    if ac != "":
                        pro = (ac, de, sq)  # 蛋白名 描述 序列
                        pl.ac.append(ac)
                        pl.de.append(de)
                        pl.sq.append(sq)
                        if using_decoy:
                            decoy_pro = CFunctionX12.GDP(ac, sq)
                            pl.ac.append(decoy_pro[0])
                            pl.de.append(decoy_pro[1])
                            pl.sq.append(decoy_pro[2])
                    line = line[1:]
                    ind = line.find(' ')
                    if ind >= 0:
                        ac = line[:ind]
                    else:
                        ac = line
                    de = line[ind + 1:]
                    sq = ""
                else:
                    sq += line
        if ac != "":
            pl.ac.append(ac)
            pl.de.append(de)
            pl.sq.append(sq)
            if using_decoy:         # 若参数self.dp.myCFG.using_decoy选择了TRUE，将反库蛋白加入protein_list中
                decoy_pro = CFunctionX12.GDP(ac, sq)
                pl.ac.append(decoy_pro[0])
                pl.de.append(decoy_pro[1])
                pl.sq.append(decoy_pro[2])
        return pl

    @staticmethod
    def GDP(ac, sq):         # 私有函数，生成反库蛋白
        new_ac = "REV_" + ac
        new_sq = sq[::-1]
        dp = (new_ac, "", new_sq)  # 蛋白名 描述 序列
        return dp

    @staticmethod
    def APr(pr:CProtein, sq:str, name:str, de:str, using_decoy=True):
        if name in pr.ac or len(name) == 0:
            CLogging.LogWarning(name)
        else:
            pr.ac.append(name)
            pr.de.append(de)
            pr.sq.append(sq)
            if using_decoy:
                # 若参数self.dp.myCFG.using_decoy选择了TRUE，将反库蛋白加入protein_list中
                decoy_pro = CFunctionX12.GDP(name, sq)
                pr.ac.append(decoy_pro[0])
                pr.de.append(decoy_pro[1])
                pr.sq.append(decoy_pro[2])

class CFunctionX13:

    @staticmethod
    def SavePkl(data, path, name):
        if not os.path.exists(path):
            os.makedirs(path)
        file_path = os.path.join(path, name)
        with open(file_path, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def LoadPkl(path, name):
        file_path = os.path.join(path, name)
        if not os.path.exists(file_path):
            CLogging.LogForLackPath(file_path)
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        return data

    @staticmethod
    def get_all_pkl_files(folder, pro_pkl_file, mod_pkl_file):
        pkl_file_list, ind_file_list = [], []
        for p in os.listdir(folder):
            if not p.endswith('.pkl'): continue
            if p == pro_pkl_file: continue
            if p == mod_pkl_file: continue
            if p.endswith("_ind.pkl"):
                #ind_file_list.append(os.path_mgf.join(folder, p))
                pass
            else:
                pkl_file_list.append(os.path.join(folder, p))
                ind_file_list.append(os.path.join(folder, p[:-4] + "_ind.pkl"))
        assert(len(pkl_file_list) == len(ind_file_list))
        return pkl_file_list, ind_file_list

class CFunctionX14:

    def __init__(self, inputDP):
        self.dp = inputDP

    def LoadSpectrum(self, path_mgf:str, path_pfind_id:str):  # 从mgf文件中提取title、mass、charge信息
        spectrum_list = []
        title_set = CFunctionX15.LoadPFindID(path_pfind_id)
        print("[Info]#PSMs in pFind is {}".format(len(title_set)))
        #title_set2 = self._load_single_spectrum(conf.single_res_path)
        #print("[Info]#PSMs in this workflow is {}".format(len(title_set2)))
        #title_set = title_set | title_set2
        title = ""
        mass = 0.0
        charge = 0
        peak_list = []
        with open(path_mgf, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line[0] == 'T':
                    title = line[6:]
                elif line[0] == 'P':
                    mass = float(line[8:])
                elif line[0] == 'C':
                    if line[-1] == '+':
                        line = line[:-1]
                    charge = int(line[7])
                elif line[0] >= '1' and line[0] <= '9':
                    ind = line.find(' ')
                    p_mz = float(line[:ind])
                    p_inten = float(line[ind+1:])
                    peak_list.append(CPeak(p_mz, p_inten))
                elif line[0] == 'E':
                    mass = mass * charge - (charge - 1) * ATOM_MASS_P
                    if len(peak_list) > self.dp.myCFG.D3_NUMBER_SELECT_PEAK:
                        cmpfun = operator.attrgetter('inten')
                        peak_list.sort(key=cmpfun, reverse=True)
                        peak_list = peak_list[:self.dp.myCFG.D3_NUMBER_SELECT_PEAK]
                    cmpfun = operator.attrgetter("mz")
                    peak_list.sort(key=cmpfun)
                    spec = CSpectrum(title, mass, charge, peak_list)
                    if len(spec.peaks) > 0 and title not in title_set:
                        spectrum_list.append(spec)
                    title = ""
                    mass = 0.0
                    charge = 0
                    peak_list = []
        return spectrum_list

    def GSMI(self, frw82frnw, G, T):

        num = int((T - G) * self.dp.myCFG.D12_MULTI_MASS) + 1
        index_list = [-1 for i in range(num)]
        for i in range(len(frw82frnw)):
            p = frw82frnw[i]
            index = int((p.mass - G) * self.dp.myCFG.D12_MULTI_MASS)
            if index >= num:
                index = num - 1
            if index_list[index] == -1:
                index_list[index] = i
        end_val = len(frw82frnw)
        for i in range(num)[::-1]:
            if index_list[i] == -1:
                index_list[i] = end_val
            else:
                end_val = index_list[i]
        return index_list

class CFunctionX15:

    @staticmethod
    def LoadPFindID(pfind_path, q_value_t=0.01):
        Q_INDEX = 4
        AC_INDEX = 12
        FLAG_INDEX = 15
        T_INDEX = 0
        title_set = set()
        if pfind_path == "" or not os.path.exists(pfind_path):
            return title_set
        with open(pfind_path) as f:
            lines = f.readlines()
        for i in range(1, len(lines)):
            tmp = lines[i].strip().split('\t')
            q_value = float(tmp[Q_INDEX])
            if q_value > q_value_t: break
            is_target = (tmp[FLAG_INDEX] == "target")
            if is_target: title_set.add(tmp[T_INDEX])
        return title_set

    @staticmethod
    def LoadSingleID(single_res_path):
        title_set = set()
        if single_res_path == "" or not os.path.exists(single_res_path):
            return title_set
        with open(single_res_path) as f:
            lines = f.readlines()
        for i in range(1, len(lines)):
            tmp = lines[i].strip().split('\t')
            title = tmp[0]
            title_set.add(title)
        return title_set

class CFunctionX16:

    @staticmethod
    def WritePLabel(psm_list, ofir382, folder, dp:CDataPack):
        pptui999_to_psm_list = {}
        for psm in psm_list:
            title = psm.spec.title
            title = title[:title.find('.')]
            if title not in pptui999_to_psm_list:
                pptui999_to_psm_list[title] = []
            pptui999_to_psm_list[title].append(psm)

        for mgf_path in ofir382:
            if mgf_path.find('/') >= 0:
                pptui999 = mgf_path[mgf_path.rfind('/')+1:]
            else:
                pptui999 = mgf_path[mgf_path.rfind('\\')+1:]
            pptui999 = pptui999[:pptui999.rfind('.')]
            if pptui999.find("_") >= 0:
                pptui999 = pptui999[:pptui999.rfind('_')]
            CFunctionX16.WriteOnepLabel(pptui999_to_psm_list[pptui999], mgf_path, folder, dp)

    @staticmethod
    def WriteOnepLabel(psm_list, mgf_path, folder, dp:CDataPack, xlink_flag="UAA"):
        if mgf_path.find('/') >= 0:
            pptui999 = mgf_path[mgf_path.rfind('/')+1:]
        else:
            pptui999 = mgf_path[mgf_path.rfind('\\')+1:]
        pptui999 = pptui999[:pptui999.rfind('.')]
        if pptui999.find("_") >= 0:
            pptui999 = pptui999[:pptui999.rfind('_')]
        path = os.path.join(folder, pptui999+".plabel")
        with open(path, 'w') as fw:
            fw.write("[FilePath]\n")
            fw.write("File_Path="+mgf_path+'\n')
            fw.write("[Modification]\n")
            mod_idx = 0
            for mod_name in dp.myINI.DIC_SET_MOD.keys():
                mod_idx += 1
                fw.write(str(mod_idx) + '=' + mod_name + '\n')
            fw.write("[xlink]\n")
            fw.write("xlink=" + xlink_flag + '\n')
            fw.write("[Total]\n")
            fw.write("total=" + str(len(psm_list)) + '\n')
            for i, psm in enumerate(psm_list):
                pep = psm.peptide
                sq1, mods1, site1 = pep.alpha_sq, pep.alpha_peptide.mods, pep.alpha_site
                sq2, mods2, site2 = pep.beta_sq, pep.beta_peptide.mods, pep.beta_site
                if type(site1) == list:
                    site1 = site1[0] + 1
                else:
                    site1 += 1
                site2 += 1
                fw.write("[Spectrum" + str(i+1) + "]\n")
                fw.write("name=" + psm.spec.title + '\n')
                fw.write("pep1=3 " + str(site1) + ' ' + str(site2) + ' ' + sq1 + ' ' + str(pep.score) + ' ' + sq2 + ' 1 ')
                mod_str = CFunctionX16.__captainGetMod2Str(mods1, 0, dp.myINI.DIC_SET_MOD)
                mod_str += CFunctionX16.__captainGetMod2Str(mods2, len(sq1) + 3, dp.myINI.DIC_SET_MOD)
                fw.write(mod_str + '\n')

    @staticmethod
    def WriteSinglePSMs(psm_list, path, is_write_candidate=False):
        fw = open(path, 'w')
        fw.write("Title\tMass\tSingle peptide sequence\tSingle peptide modification mass\tSingle peptide site\t" +
                 "Score\tDelta score\tRank\tTarget-Decoy\n")
        for p in psm_list:
            title = p.spec.title
            mass = p.spec.mass
            delta_score = p.peptide.score
            if len(p.peptide_list) > 1:
                delta_score -= p.peptide_list[1].score
            is_target = (p.peptide.pro_index == 0)
            fw.write(title + '\t' + str(mass) + '\t' + p.peptide.peptide_sq + '\t' + str(p.peptide.open_mass) +
                     '\t' + str(p.peptide.open_site + 1) + '\t' + str(p.peptide.score) + '\t' + str(delta_score) + '\t1\t' +
                     str(is_target) + '\n')
        fw.close()

    @staticmethod
    def WritePSMs(psm_list, path, conf:CConfig, is_write_candidate=False):
        protein_list = CFunctionX13.LoadPkl(conf.A4_PATH_FASTA_EXPORT, NAME_PROTEIN_PKL)
        fw = open(path, 'w')
        fw.write("Title\tMass\tAlpha peptide sequence\tAlpha peptide modification\tAlpha site\tAlpha mass\t" +
                "Beta peptide sequence\tBeta peptide modification\tBeta site\tBeta mass\t" +
                "Match tolerance (ppm/Da)\tProteins\tScore\tDelta score\tRank\t" +
                 "Alpha score\tAlpha mean fragment mass tol\t" +
                 "Alpha std fragment mass tol\tAlpha ion ratio\tAlpha inten ratio\t" +
                 "Alpha b-score\tAlpha y-score\tAlpha b/y-score\tAlpha sq len\t" +
                 "Beta score\tBeta mean fragment mass tol\t" +
                 "Beta std fragment mass tol\tBeta ion ratio\tBeta inten ratio\t" +
                 "Beta b-score\tBeta y-score\tBeta b/y-score\tBeta sq len\n")
        for p in psm_list:
            title = p.spec.title
            mass = p.spec.mass
            delta_score = p.peptide.alpha_score.match_score + p.peptide.beta_score.match_score
            if len(p.peptide_list) > 1:
                delta_score -= (p.peptide_list[1].alpha_score.match_score + p.peptide_list[1].beta_score.match_score)
            delta_score = abs(delta_score)
            CFunctionX16.__captainWriteOnePSM(title, mass, p.peptide, conf.C1_TYPE_TOL_PRECURSOR, fw, protein_list, delta_score)
            if is_write_candidate:
                for i in range(1,len(p.peptide_list)):
                    CFunctionX16.__captainWriteOnePSM(title, mass, p.peptide_list[i], conf.C1_TYPE_TOL_PRECURSOR, fw, protein_list, 0.0, i + 1)
        fw.close()

    @staticmethod
    def __captainGetMod2Str(mods, start_site, mod_2_index):
        mod_str = ""
        for m in mods:
            site = m.site + start_site
            name = m.mod_name
            index = list(mod_2_index.keys()).index(name)+1
            if name.find('N-term') >= 0: site -= 1
            elif name.find('C-term') >= 0: site += 1
            mod_str += str(site) + ',' + str(index) + ' '
        return mod_str

    @staticmethod
    def __captainGetModStr(mods):
        mod_str = ""
        for m in mods:
            mod_str += str(m.site) + "@" + str(m.mod_name) + ";"
        if len(mod_str) > 0 and mod_str[-1] == '\t':
            mod_str = mod_str[:-1]
        return mod_str

    @staticmethod
    def __captainWriteOnePSM(title, mass, pep:CLinkPeptide, is_ppm_pre, fw, protein_list:CProtein, delta_score = 0.0, rank=1):
        sq1, mods1, site1, mass1, score_f1 = pep.alpha_sq, pep.alpha_peptide.mods, pep.alpha_site, pep.alpha_peptide.mass, pep.alpha_score
        sq2, mods2, site2, mass2, score_f2 = pep.beta_sq, pep.beta_peptide.mods, pep.beta_site, pep.beta_peptide.mass, pep.beta_score
        mods1 = CFunctionX16.__captainGetModStr(mods1)
        mods2 = CFunctionX16.__captainGetModStr(mods2)
        if type(site1) == list:
            if len(site1) > 1:
                new_site1 = copy.deepcopy(site1)
                for i in range(len(site1)):
                    new_site1[i] = str(site1[i] + 1)
                site1 = ";".join(new_site1)
            else:
                site1 = site1[0] + 1
        else:
            site1 += 1# start from 1 not 0
        site2 += 1
        ac_list = [protein_list.ac[i] for i in pep.alpha_peptide.pro_index_list]
        ac_str = ";".join(ac_list)
        if is_ppm_pre: match_tolerance = (mass - mass1 - mass2 - ATOM_MASS_P) * 1e6 / mass
        else: match_tolerance = (mass - mass1 - mass2 - ATOM_MASS_P)
        fw.write(title + '\t' + "%.6f"%(mass) + '\t' + sq1 + '\t' + mods1 + '\t' + str(site1) + '\t' + "%.6f"%(mass1) + '\t' +
                 sq2 + '\t' + mods2 + '\t' + str(site2) + '\t' + "%.6f"%(mass2) + '\t' +
                 "%.6f"%(match_tolerance) + '\t' + ac_str + '\t' + "%.6f"%(pep.alpha_score.score) + '\t' + "%.6f"%(delta_score) + '\t' + str(rank) + '\t' + score_f1.get_string() + '\t' + str(len(sq1)) + '\t' + score_f2.get_string() + '\t' + str(len(sq2)) + '\n')

