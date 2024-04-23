from MSSystem import MODEL_PATH
from MSDataC import CDataPack

from MSFunctionD import CFunctionX6

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
import math

class CFunctionX23:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def Rerank(self, psm_list, result_file_list, feature_write_path):

        td_list, record_filter_num = CFunctionX6.updatePSMFDRvalue(psm_list)

        ori_psm_list = copy.deepcopy(psm_list)

        tmp_filter_num = record_filter_num
        for k in range(self.dp.RerankTime):

            self.__caEIORPQJPFDSjfdsoijoiewq(psm_list, feature_write_path + '.test')

            td_list, tmp_filter_num = CFunctionX6.updatePSMFDRvalue(psm_list)
            print(tmp_filter_num)
            target_list = []
            decoy_list = []
            for i in range(len(td_list)):
                if not psm_list[i].decoy_flag:
                    if td_list[i] < self.dp.myCFG.E1_FDR_PSM:
                        target_list.append(psm_list[i])
                else:
                    decoy_list.append(psm_list[i])
            if (len(target_list) == 0) or (len(decoy_list) == 0):
                break
            # if k == 0:
            #     self._draw_feature_density_plot(target_list, decoy_list, round=k)
            self.__caEIORPQJPFDSjfdsoijoiewq(target_list + decoy_list, feature_write_path + '.train')

            self.__captainGetDataScale(feature_write_path + '.test', feature_write_path + '.test' + '.scale')
            self.__captainGetDataScale(feature_write_path + '.train', feature_write_path + '.train' + '.scale')

            score_write_path = feature_write_path + '.score'

            self.__captainSVMTrain(feature_write_path + '.test' + '.scale', feature_write_path + '.train' + '.scale', score_write_path)

            score_list = self.__captainGetSVMScore(score_write_path)
            score_list = [float(score) for score in score_list]

            assert(len(score_list) == len(psm_list))

            self.__captainUpdatePSMOrder(psm_list, score_list)

            td_list, check_filter_num = CFunctionX6.updatePSMFDRvalue(psm_list)
            if check_filter_num == tmp_filter_num:
                break

        td_list, final_filter_num = CFunctionX6.updatePSMFDRvalue(psm_list)
        if final_filter_num > record_filter_num:
            psm_dict = {}
            for psm in psm_list:
                psm_dict[psm.title] = psm.link_peptide.svm_score
            self.__captainUpdatePSMScore(result_file_list, psm_dict)
            return True
        else:
            psm_dict = {}
            for psm in ori_psm_list:
                psm_dict[psm.title] = psm.link_peptide.svm_score
            self.__captainUpdatePSMScore(result_file_list, psm_dict)
            return False

    def __captaFJEOIWJFEOIWASD(self, target_list, decoy_list, round):
        total_score_t = np.array([(target[1].alpha_score.score + target[1].beta_score.score) for target in target_list])
        total_score_d = np.array([(decoy[1].alpha_score.score + decoy[1].beta_score.score) for decoy in decoy_list])
        plt.figure(figsize=(10,10))
        plt.title("Total score")
        plt.hist(total_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(total_score_t, color="r")
        plt.hist(total_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(total_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        os.path.join(self.dp.myCFG.A5_PATH_RESULT_EXPORT)
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\total_score-{round}.png")

        delta_score_t = np.array([(target[-2]) for target in target_list])
        delta_score_d = np.array([(decoy[-2]) for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("Delta score")
        plt.hist(delta_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(delta_score_t, color="r")
        plt.hist(delta_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(delta_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\delta_score-{round}.png")

        alpha_score_t = np.array([target[1].alpha_score.score for target in target_list])
        alpha_score_d = np.array([decoy[1].alpha_score.score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha score")
        plt.hist(alpha_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(alpha_score_t, color="r")
        plt.hist(alpha_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(alpha_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_score-{round}.png")

        alpha_con_score_t = np.array([target[1].alpha_score.continue_score for target in target_list])
        alpha_con_score_d = np.array([decoy[1].alpha_score.continue_score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha continuous score")
        plt.hist(alpha_con_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(alpha_con_score_t, color="r")
        plt.hist(alpha_con_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(alpha_con_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_continuous_score-{round}.png")

        alpha_con_y_score_t = np.array([target[1].alpha_score.continue_y_score for target in target_list])
        alpha_con_y_score_d = np.array([decoy[1].alpha_score.continue_y_score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha y ions continuous score")
        plt.hist(alpha_con_y_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(alpha_con_y_score_t, color="r")
        plt.hist(alpha_con_y_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(alpha_con_y_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_continuous_y_score-{round}.png")

        alpha_con_b_score_t = np.array([target[1].alpha_score.continue_b_score for target in target_list])
        alpha_con_b_score_d = np.array([decoy[1].alpha_score.continue_b_score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha b ions continuous score")
        plt.hist(alpha_con_b_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(alpha_con_b_score_t, color="r")
        plt.hist(alpha_con_b_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(alpha_con_b_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_continuous_b_score-{round}.png")

        beta_score_t = np.array([target[1].beta_score.score for target in target_list])
        beta_score_d = np.array([decoy[1].beta_score.score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("beta score")
        plt.hist(beta_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(beta_score_t, color="r")
        plt.hist(beta_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(beta_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\beta_score-{round}.png")

        beta_con_score_t = np.array([target[1].beta_score.continue_score for target in target_list])
        beta_con_score_d = np.array([decoy[1].beta_score.continue_score for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("beta continuous score")
        plt.hist(beta_con_score_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(beta_con_score_t, color="r")
        plt.hist(beta_con_score_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(beta_con_score_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\beta_continuous_score-{round}.png")

        len_t = np.array([len(target[1].alpha_peptide.sq) for target in target_list])
        len_d = np.array([len(decoy[1].alpha_peptide.sq) for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide length")
        plt.hist(len_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(len_t, color="r")
        plt.hist(len_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(len_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_pep_length-{round}.png")

        log_len_t = np.array([math.log(min(len(target[1].alpha_peptide.sq),len(target[1].beta_peptide.sq))) for target in target_list])
        log_len_d = np.array([math.log(min(len(decoy[1].alpha_peptide.sq),len(decoy[1].beta_peptide.sq))) for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("log minimum peptide length")
        plt.hist(log_len_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(log_len_t, color="r")
        plt.hist(log_len_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")  # 绘制直方图
        sns.kdeplot(log_len_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\log_minimum_length-{round}.png")

        inten_t = np.array([target[1].alpha_score.inten_ratio for target in target_list])
        inten_d = np.array([decoy[1].alpha_score.inten_ratio for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide match intensity")
        plt.hist(inten_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(inten_t, color="r")
        plt.hist(inten_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(inten_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_peptide_match_intensity-{round}.png")

        ion_t = np.array([target[1].alpha_score.ion_ratio for target in target_list])
        ion_d = np.array([decoy[1].alpha_score.ion_ratio for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide match ions")
        plt.hist(ion_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(ion_t, color="r")
        plt.hist(ion_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(ion_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_peptide_match_ions-{round}.png")

        mean_error_t = np.array([target[1].alpha_score.mean_fra_me for target in target_list])
        mean_error_d = np.array([decoy[1].alpha_score.mean_fra_me for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide match error mean")
        plt.hist(mean_error_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(mean_error_t, color="r")
        plt.hist(mean_error_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(mean_error_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_peptide_match_error_mean-{round}.png")

        std_error_t = np.array([target[1].alpha_score.std_fra_me for target in target_list])
        std_error_d = np.array([decoy[1].alpha_score.std_fra_me for decoy in decoy_list])
        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide match error standard deviation")
        plt.hist(std_error_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(std_error_t, color="r")
        plt.hist(std_error_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(std_error_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_peptide_match_error_std-{round}.png")

        miss_clv_t = []
        for i, (title, link_pep, decoy_flag, ranker, delta_score, svm_score) in enumerate(target_list):
            enzyme_rule = self.dp.myCFG.B1_NAME_ENZYME.split(" ")[1:3]
            alpha_miss_clv_num = 0
            if enzyme_rule[1] == "C":
                for j in range(len(link_pep.alpha_peptide.sq) - 1):
                    if link_pep.alpha_peptide.sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1
            else:
                for j in range(1, len(link_pep.alpha_peptide.sq)):
                    if link_pep.alpha_peptide.sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1
            miss_clv_t.append(alpha_miss_clv_num/math.log(len(link_pep.alpha_peptide.sq)))
        miss_clv_d = []
        for i, (title, link_pep, decoy_flag, ranker, delta_score, svm_score) in enumerate(decoy_list):
            enzyme_rule = self.dp.myCFG.B1_NAME_ENZYME.split(" ")[1:3]
            alpha_miss_clv_num = 0
            if enzyme_rule[1] == "C":
                for j in range(len(link_pep.alpha_peptide.sq) - 1):
                    if link_pep.alpha_peptide.sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1
            else:
                for j in range(1, len(link_pep.alpha_peptide.sq)):
                    if link_pep.alpha_peptide.sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1
            miss_clv_d.append(alpha_miss_clv_num/math.log(len(link_pep.alpha_peptide.sq)))

        plt.figure(figsize=(10, 10))
        plt.title("alpha peptide length de miss clv num")
        plt.hist(miss_clv_t, bins=30, density=True, alpha=0.5, color="r", label="target")
        sns.kdeplot(miss_clv_t, color="r")
        plt.hist(miss_clv_d, bins=30, density=True, alpha=0.5, color="b", label="decoy")
        sns.kdeplot(miss_clv_d, color="b")
        plt.legend(['Target', 'Target', 'Decoy', 'Decoy'])
        plt.savefig(f"E:\\dandan\\eFSY\\TRX-10\\REP3\\result_efsy(pf)\\td_image\\alpha_peptide_clv-{round}.png")

    def __caEIORPQJPFDSjfdsoijoiewq(self, psm_list, write_path):
        fw = open(write_path, 'w')
        for i, one_psm in enumerate(psm_list):
            alpha_score, beta_score = one_psm.link_peptide.alpha_score, one_psm.link_peptide.beta_score
            len1 = len(one_psm.link_peptide.alpha_peptide.sq)
            len2 = len(one_psm.link_peptide.beta_peptide.sq)
            total_score = alpha_score.score
            total_score += beta_score.score
            class_flag = 1
            if one_psm.decoy_flag: class_flag = -1
            alpha_sq = one_psm.link_peptide.alpha_peptide.sq
            enzyme_rule = self.dp.myCFG.B1_NAME_ENZYME.split(" ")[1:3]
            alpha_miss_clv_num = 0
            if enzyme_rule[1] == "C":
                for j in range(len(alpha_sq)-1):
                    if alpha_sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1
            else:
                for j in range(1, len(alpha_sq)):
                    if alpha_sq[j] in enzyme_rule[0]:
                        alpha_miss_clv_num += 1

            fw.write(f"{class_flag} 1:{total_score} 2:{one_psm.delta_score} 3:{alpha_score.score} "
                     f"4:{alpha_score.continue_score} 5:{alpha_score.inten_ratio} 6:{alpha_score.ion_ratio} "
                     f"7:{alpha_miss_clv_num / math.log(len1)} 8:{one_psm.link_peptide.svm_score}\n")

        fw.close()

    def __captainGetDataScale(self, test_file, test_scale_file):

        cmd_code = 'cd ' + os.getcwd() + "\\libsvm\\windows" + ' && '
        test_scale_code = cmd_code + "svm-scale.exe" + ' ' + test_file + ' > ' + test_scale_file

        os.system(test_scale_code)

    def __captainSVMTrain(self, test_scale_file, train_scale_file, score_write_path):

            model_file = train_scale_file + '.model'

            cmd_code = 'cd ' + os.getcwd() + '\\' + "\\libsvm\\windows" + ' && '
            train_code = cmd_code + "svm-train.exe" + ' -b 1 -h 0  -c ' + str(1) + ' ' + train_scale_file + ' ' + model_file
            predict_code = cmd_code + "svm-predict.exe" + ' -b 1 ' + test_scale_file + ' ' + model_file + ' ' + score_write_path

            os.system(train_code)
            os.system(predict_code)

    def __captainSVMClassify(self, test_path, model_path, output_path):
        os.system(os.getcwd() + "\\libsvm\\windows\\svm-predict.exe -b 1 " + test_path + " " + model_path + " " + output_path)

    def __captainUpdatePSMScore(self, result_file_list, psm_dict):
        for file_path in result_file_list:
            write_path = file_path + "_new.txt"
            tmp_list = []
            fw = open(write_path, 'w')
            head_line = ""
            with open(file_path, 'r') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    if line_num <= 1:
                        head_line = line.strip().split('\t')
                        head_line.insert(13, "Svm score")
                        new_head_line = "\t".join(head_line)
                        continue
                    tmp = line.strip().split('\t')
                    tmp.insert(13, str(psm_dict[tmp[0]]))
                    tmp_list.append((psm_dict[tmp[0]], "\t".join(tmp)))
            tmp_list.sort(key=lambda tup: tup[0], reverse=True)
            fw.write(new_head_line + '\n')
            for (_, line) in tmp_list:
                fw.write(line + '\n')
            fw.close()

    def __captainUpdatePSMOrder(self, psm_list, score_list):

        for i, one_psm in enumerate(psm_list):
            psm_list[i].link_peptide.svm_score = score_list[i]

        psm_list.sort(key=lambda x: x.link_peptide.svm_score, reverse=True)

    def __captainGetSVMScore(self, score_path):
        score_list = []
        pos_flag = 1
        with open(score_path) as f:
            for line in f:
                line = line.strip()
                if line == "": continue
                if line.startswith("labels"):
                    tmp = line.split(' ')
                    if tmp[2] == '1': pos_flag = 2
                    continue
                tmp = line.split(' ')
                score_list.append(str(tmp[pos_flag]))
        return score_list