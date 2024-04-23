# from MSData import CModSite, CSinglePeptide, CPeptide, CLinkPeptide, CSingle_Score_Feature
from MSSystem import XLINK_PEP, XLINK_PSM, SINGLE_PEP, SINGLE_PSM
from MSDataC import CDataPack
from MSDataPeptide import CModSite
from MSDataA import CPeptideIDFeature, CPeptideIDInfo, CLinkPeptideIDInfo, CCrosslinkPSMID, CSinglePeptideIDInfo, CSinglePSMID


import time
import os


class CFunctionX5:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def GetTaskId(self, need_time=True):  # 给每次进行重打分的结果文件赋一个特殊的值作为文件名，防止重复覆盖
        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        if len(ofir382) > 0:
            pptui999 = ofir382[0]
        else:
            pptui999 = self.dp.myCFG.A2_PATH_MS2
        if pptui999.find('/') >= 0:
            pptui999 = pptui999.split('/')[-1]
        else:
            pptui999 = pptui999.split('\\')[-1]
        pptui999 = pptui999[:pptui999.rfind('.')]
        if need_time:
            taskid = pptui999 + "-" + str(int(time.time()))
        else:
            taskid = pptui999
        return taskid

    def FindSingleResFile(self):
        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        iitu382 = set()
        for file in ofir382:
            if file.find('/') >= 0:
                file = file.split('/')[-1][:-4]
            else:
                file = file.split('\\')[-1][:-4]
            iitu382.add(file + "_singlePeptide.txt")
        folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT
        file_list = []
        for file in os.listdir(folder):
            if file in iitu382: file_list.append(os.path.join(folder, file))
        return file_list

    def FindResFile(self, cand_flag = False, label_flag = False, rerank_flag = False):
        ofir382 = self.dp.myCFG.A2_PATH_MS2.split(";")
        iitu382 = set()
        for file in ofir382:
            if file.find('/') >= 0:
                file = file.split('/')[-1][:-4]
            else:
                file = file.split('\\')[-1][:-4]
            if cand_flag:
                file += "_cand"
            if not cand_flag and label_flag:
                index = file.rfind('_')
                if index >= 0:
                    file = file[:index]
            if not label_flag:
                if rerank_flag:
                    file += ".txt_new.txt"
                else:
                    file += ".txt"
            else:
                file += ".plabel"
            iitu382.add(file)
        folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT
        file_list = []
        for file in os.listdir(folder):
            if not label_flag and not rerank_flag and not file.endswith(".txt"):
                continue
            if label_flag and not file.endswith(".plabel"):
                continue
            if rerank_flag and not file.endswith(".txt_new.txt"):
                continue
            if file in iitu382:
                file_list.append(os.path.join(folder, file))
        return file_list

    def __captainGetModS(self, mod_str):
        mod_list = []
        for m in mod_str.split(';'):
            if m.find('@') < 0: continue
            site, name = m.split('@')
            site = int(site)
            one_mod = CModSite(name, site, self.dp.myINI.DIC_MOD[name].mass)
            mod_list.append(one_mod)
        return mod_list

    def _get_mod2mass(self, mod2mod):
        mod2mass = {}
        for name in mod2mod:
            mod2mass[name] = mod2mod[name].mass
        return mod2mass

    def LoadOpenUAASinglePSM(self, path):
        psm_list = []
        with open(path, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num <= 1: continue
                title, mass, sq, mod_mass, site, score, delta_score, rank, is_target = line.strip().split('\t')
                mass = float(mass)
                mod_mass = float(mass)
                site = int(site)
                score = float(score)
                delta_score = float(delta_score)
                rank = int(rank)
                decoy_flag = (is_target != "True")
                pro_index = 0
                if decoy_flag: pro_index = -1
                peptide = CSinglePeptideIDInfo(sq, mass, mod_mass, pro_index, score)
                psm_list.append(CSinglePSMID(title, peptide, decoy_flag, rank, delta_score))
        return psm_list

    def LoadOpenUAACrosslinkPSM(self, path):
        psm_list = []
        with open(path, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num <= 1: continue
                tmp = line.strip().split('\t')
                score = float(tmp[12])
                delta_score = float(tmp[13])
                title = tmp[0]
                sq1 = tmp[2]
                mods1 = tmp[3]
                if ";" in tmp[4]:
                    site1 = int(tmp[4].split(";")[0]) - 1
                else:
                    site1 = int(tmp[4]) - 1
                mass1 = float(tmp[5])
                sq2 = tmp[6]
                mods2 = tmp[7]
                site2 = int(tmp[8]) - 1
                mass2 = float(tmp[9])
                ac_list = tmp[11]
                decoy_flag = True
                for ac in ac_list.split(';'):
                    if not ac.startswith("REV_"):
                        decoy_flag = False
                        break
                ranker = int(tmp[14])
                score1 = float(tmp[15])
                mean_fra_me1 = float(tmp[16])
                std_fra_me1 = float(tmp[17])
                ion_ratio1 = float(tmp[18])
                inten_ratio1 = float(tmp[19])
                continue_b_score1 = float(tmp[20])
                continue_y_score1 = float(tmp[21])
                continue_score1 = float(tmp[22])
                len1 = int(tmp[23])
                score2 = float(tmp[24])
                mean_fra_me2 = float(tmp[25])
                std_fra_me2 = float(tmp[26])
                ion_ratio2 = float(tmp[27])
                inten_ratio2 = float(tmp[28])
                continue_b_score2 = float(tmp[29])
                continue_y_score2 = float(tmp[30])
                continue_score2 = float(tmp[31])
                len2 = int(tmp[32])
                alpha_score = CPeptideIDFeature(score1, mean_fra_me1, std_fra_me1, ion_ratio1, inten_ratio1, continue_b_score1, continue_y_score1, continue_score1)
                beta_score = CPeptideIDFeature(score2, mean_fra_me2, std_fra_me2, ion_ratio2, inten_ratio2, continue_b_score2, continue_y_score2, continue_score2)
                pep1 = CPeptideIDInfo(sq1, mass1, self.__captainGetModS(mods1), [])
                pep2 = CPeptideIDInfo(sq2, mass2, self.__captainGetModS(mods2), [])
                link_peptide = CLinkPeptideIDInfo(pep1, pep2, site1, site2, score, alpha_score, beta_score)
                psm_list.append(CCrosslinkPSMID(title, link_peptide, decoy_flag, ranker, delta_score))
        return psm_list

    def LoadOpenUAARerankCrosslinkPSM(self, path):
        psm_list = []
        with open(path, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num <= 1: continue
                tmp = line.strip().split('\t')
                score = float(tmp[12])
                svm_score = float(tmp[13])
                delta_score = float(tmp[14])
                title = tmp[0]
                sq1 = tmp[2]
                mods1 = tmp[3]
                if ";" in tmp[4]:
                    site1 = int(tmp[4].split(";")[0]) - 1
                else:
                    site1 = int(tmp[4]) - 1
                mass1 = float(tmp[5])
                sq2 = tmp[6]
                mods2 = tmp[7]
                site2 = int(tmp[8]) - 1
                mass2 = float(tmp[9])
                ac_list = tmp[11]
                decoy_flag = True
                for ac in ac_list.split(';'):
                    if not ac.startswith("REV_"):
                        decoy_flag = False
                        break
                ranker = int(tmp[15])
                score1 = float(tmp[16])
                mean_fra_me1 = float(tmp[17])
                std_fra_me1 = float(tmp[18])
                ion_ratio1 = float(tmp[19])
                inten_ratio1 = float(tmp[20])
                continue_b_score1 = float(tmp[21])
                continue_y_score1 = float(tmp[22])
                continue_score1 = float(tmp[23])
                len1 = int(tmp[24])
                score2 = float(tmp[25])
                mean_fra_me2 = float(tmp[26])
                std_fra_me2 = float(tmp[27])
                ion_ratio2 = float(tmp[28])
                inten_ratio2 = float(tmp[29])
                continue_b_score2 = float(tmp[30])
                continue_y_score2 = float(tmp[31])
                continue_score2 = float(tmp[32])
                len2 = int(tmp[33])
                alpha_score = CPeptideIDFeature(score1, mean_fra_me1, std_fra_me1, ion_ratio1, inten_ratio1, continue_b_score1, continue_y_score1, continue_score1)
                beta_score = CPeptideIDFeature(score2, mean_fra_me2, std_fra_me2, ion_ratio2, inten_ratio2, continue_b_score2, continue_y_score2, continue_score2)
                pep1 = CPeptideIDInfo(sq1, mass1, self.__captainGetModS(mods1), [])
                pep2 = CPeptideIDInfo(sq2, mass2, self.__captainGetModS(mods2), [])
                link_peptide = CLinkPeptideIDInfo(pep1, pep2, site1, site2, score, alpha_score, beta_score)
                link_peptide.svm_score = svm_score
                psm_list.append(CCrosslinkPSMID(title, link_peptide, decoy_flag, ranker, delta_score))
        return psm_list


class CFunctionX6:

    def __init__(self, inputDP):

        self.dp = inputDP

    def calculateFDRvalueBySVM(self, psm_list):

        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        for (title, pep, decoy, ranker) in psm_list:
            if decoy:
                decoy_num += 1
            else:
                target_num += 1
        print("[Info] #PSMs is {0}, #Target is {1} and #Decoy is {2}".format(all_num, target_num, decoy_num))
        cur_target_num, cur_decoy_num = 0, 0
        if len(psm_list) == 0: return 0.0
        score_t = psm_list[0][1].svm_score + 1.0
        if target_num > 0 and decoy_num * 1.0 / target_num <= self.dp.myCFG.E1_FDR_PSM:
            score_t = psm_list[-1][1].svm_score
            return score_t
        for (title, pep, decoy, ranker) in psm_list[::-1]:
            if decoy:
                cur_decoy_num += 1
            else:
                cur_target_num += 1
            if target_num == cur_target_num:
                fdr_value = 100
            else:
                fdr_value = (decoy_num - cur_decoy_num) * 1.0 / (target_num - cur_target_num)
            # print(decoy_num,cur_decoy_num, target_num, cur_target_num, fdr_value)
            if fdr_value <= self.dp.myCFG.E1_FDR_PSM:
                score_t = pep.svm_score
                return score_t
        return score_t

    def calculateFDRvalueByORI(self, psm_list):

        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        for (title, pep, decoy, ranker) in psm_list:
            if decoy:
                decoy_num += 1
            else:
                target_num += 1
        print("[Info] #PSMs is {0}, #Target is {1} and #Decoy is {2}".format(all_num, target_num, decoy_num))
        cur_target_num, cur_decoy_num = 0, 0
        if len(psm_list) == 0: return 0.0
        score_t = psm_list[0][1].score + 1.0
        if target_num > 0 and decoy_num * 1.0 / target_num <= self.dp.myCFG.E1_FDR_PSM:
            score_t = psm_list[-1][1].score
            return score_t
        for (title, pep, decoy, ranker) in psm_list[::-1]:
            if decoy:
                cur_decoy_num += 1
            else:
                cur_target_num += 1
            if target_num == cur_target_num:
                fdr_value = 100
            else:
                fdr_value = (decoy_num - cur_decoy_num) * 1.0 / (target_num - cur_target_num)
            # print(decoy_num,cur_decoy_num, target_num, cur_target_num, fdr_value)
            if fdr_value <= self.dp.myCFG.E1_FDR_PSM:
                score_t = pep.score
                return score_t
        return score_t

    @staticmethod
    def updatePSMFDRvalue(psm_list):

        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        td_list = []
        for psm in psm_list:
            if psm.decoy_flag:
                decoy_num += 1
            else:
                target_num += 1
            if target_num != 0:
                td_list.append(decoy_num / target_num)
            else:
                td_list.append(1)

        tmp_fdr = td_list[-1]
        for i in reversed(range(len(td_list))):
            if td_list[i] < tmp_fdr:
                tmp_fdr = td_list[i]
            else:
                td_list[i] = tmp_fdr

        record_filter_num = 0
        for i in range(len(td_list)):
            if td_list[i] > 0.05:
                record_filter_num = i
                break

        return td_list, record_filter_num


    def write_FDR_result(self, res_file_list, res_plabel_file_list, fdr_score, taskid, svm_status):
        if svm_status:
            check_index = 13
        else:
            check_index = 12
        pep_score, pep_line = {}, {} # 记录肽段水平的最佳匹配分数和信息
        t_p_line = {}
        head_line = ""
        title_score, title_line, title_plabel_line = {}, {}, {} # 记录谱图水平的匹配分数和信息
        pep_plabel_line = {}
        for file in res_file_list:
            with open(file) as f:
                lines = f.readlines()
            head_line = lines[0].strip()
            for i in range(1, len(lines)):
                line = lines[i].strip()
                tmp = line.split('\t')
                title = tmp[0]
                sq_key = tmp[2] + '@' + tmp[3] + '@' + tmp[4] + '@' + tmp[6]
                ac_list = tmp[11].split(';')
                is_target = False
                for ac in ac_list:
                    if not ac.startswith("REV_"):
                        is_target = True
                        break
                if not is_target: continue
                s = float(tmp[check_index])
                if s < fdr_score: break
                title_score[title] = s
                title_line[title] = line
                if sq_key not in pep_score or pep_score[sq_key] < s:
                    pep_score[sq_key] = s
                    pep_line[sq_key] = line

        for sq in pep_line:
            line = pep_line[sq]
            tmp = line.strip().split('\t')
            title = tmp[0]
            t_p_line[title] = 1

        modification_lines = []
        for file in res_plabel_file_list:
            with open(file) as f:
                lines = f.readlines()
            for i in range(len(lines)):
                if lines[i].startswith("name="):
                    title = lines[i].strip().split('=')[1].strip()
                    if title in title_line:
                        line = lines[i + 1].strip()
                        title_plabel_line[title] = line
                    if title in t_p_line:
                        line = lines[i + 1].strip()
                        pep_plabel_line[title] = line
                elif lines[i].startswith("[Modification]"):
                    modification_lines = []
                    for j in range(i + 1, len(lines)):
                        ind = lines[j].find('=')
                        if ind < 0: break
                        modification_lines.append(lines[j].strip())

        sorted_res = sorted(title_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, XLINK_PSM + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = title_line[r]
                f.write(one_line + '\n')
        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, XLINK_PSM + '-' + taskid + '.plabel')
        with open(write_path, 'w') as f:
            f.write("[FilePath]\n")
            f.write("File_Path=" + ";".join(self.dp.myCFG.A2_PATH_MS2.split(";")) + "\n")
            f.write("[Modification]\n")
            for mod_line in modification_lines:
                f.write(mod_line + '\n')
            f.write("[xlink]\n")
            f.write("xlink=UAA\n")
            f.write("[Total]\n")
            f.write("total=" + str(len(sorted_res)) + '\n')
            spec_id = 0
            for (r, s) in sorted_res:
                spec_id += 1
                one_line = title_plabel_line[r]
                f.write("[Spectrum" + str(spec_id) + "]\n")
                f.write("name=" + r + "\n")
                f.write(one_line + "\n")

        sorted_res = sorted(pep_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, XLINK_PEP + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = pep_line[r]
                f.write(one_line + '\n')

        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, XLINK_PEP + '-' + taskid + '.plabel')
        with open(write_path, 'w') as f:
            f.write("[FilePath]\n")
            f.write("File_Path=" + ";".join(self.dp.myCFG.A2_PATH_MS2.split(";")) + "\n")
            f.write("[Modification]\n")
            for mod_line in modification_lines:
                f.write(mod_line + '\n')
            f.write("[xlink]\n")
            f.write("xlink=UAA\n")
            f.write("[Total]\n")
            f.write("total=" + str(len(sorted_res)) + '\n')
            spec_id = 0
            for (r, s) in sorted_res:
                spec_id += 1
                one_line = pep_line[r]
                tmp = one_line.strip().split('\t')
                title = tmp[0]
                f.write("[Spectrum" + str(spec_id) + "]\n")
                f.write("name=" + title + "\n")
                one_line = pep_plabel_line[title]
                f.write(one_line + "\n")


class CFunctionX7:

    def __init__(self, inputDP):
        self.dp = inputDP

    def CalculateFDR(self, psm_list):
        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        for psm in psm_list:
            if psm.decoy_flag:
                decoy_num += 1
            else:
                target_num += 1
        print("[Info] #PSMs is {0}, #Target is {1} and #Decoy is {2}".format(all_num, target_num, decoy_num))
        cur_target_num, cur_decoy_num = 0, 0
        if len(psm_list) == 0: return 0.0
        score_t = psm_list[0].peptide.score + 1.0
        if target_num > 0 and decoy_num * 1.0 / target_num <= self.dp.myCFG.E1_FDR_PSM:
            score_t = psm_list[-1].peptide.score
            return score_t
        for (title, pep, decoy, ranker, delta_score) in psm_list[::-1]:
            if decoy:
                cur_decoy_num += 1
            else:
                cur_target_num += 1
            if target_num == cur_target_num:
                fdr_value = 100
            else:
                fdr_value = (decoy_num - cur_decoy_num) * 1.0 / (target_num - cur_target_num)
            # print(decoy_num,cur_decoy_num, target_num, cur_target_num, fdr_value)
            if fdr_value <= self.dp.myCFG.E1_FDR_PSM:
                score_t = pep.score
                return score_t
        return score_t

    def write_FDR_result(self, res_file_list, fdr_score, taskid):
        pep_score, pep_line = {}, {}
        head_line = ""
        title_score, title_line = {}, {}
        for file in res_file_list:
            with open(file) as f:
                lines = f.readlines()
            head_line = lines[0].strip()
            for i in range(1, len(lines)):
                line = lines[i].strip()
                tmp = line.split('\t')
                title = tmp[0]
                sq_key = tmp[2] + '@' + tmp[4] + '@' + str(int(float(tmp[3]) + 0.5))
                is_target = (tmp[8] == "True")
                if not is_target: continue
                s = float(tmp[5])
                if s <= fdr_score: break
                title_score[title] = s
                title_line[title] = line
                if sq_key not in pep_score or pep_score[sq_key] < s:
                    pep_score[sq_key] = s
                    pep_line[sq_key] = line

        sorted_res = sorted(title_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, SINGLE_PSM + '-' + taskid + '.txt')
        SINGLE_RES_PATH = write_path
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = title_line[r]
                f.write(one_line + '\n')

        sorted_res = sorted(pep_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT, SINGLE_PEP + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = pep_line[r]
                f.write(one_line + '\n')