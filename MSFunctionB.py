
from MSLogging import CLogging
from MSDataC import CDataPack
from MSDataPeptide import CPeptideInfo
from MSSystem import NAME_PROTEIN_PKL, NAME_MOD_PKL
from MSTool import tool_get_start_T

import os
import pickle

class CFunctionX1:

    def __init__(self, inputDP:CDataPack):

        self.dp = inputDP

    def gedpi(self):
        pkl_num = int((self.dp.myCFG.D7_MASS_PEP_UP - self.dp.myCFG.D6_MASS_PEP_LOW) / self.dp.myCFG.D10_INDEX_SPLIT_MASS) + 1
        hjs_list = []
        for i in range(pkl_num):
            G = self.dp.myCFG.D6_MASS_PEP_LOW + i * self.dp.myCFG.D10_INDEX_SPLIT_MASS
            T = self.dp.myCFG.D6_MASS_PEP_LOW + (i + 1) * self.dp.myCFG.D10_INDEX_SPLIT_MASS
            if T > self.dp.myCFG.D7_MASS_PEP_UP:
                T = self.dp.myCFG.D7_MASS_PEP_UP
            if G == T:
                continue
            folder = self.dp.myCFG.A4_PATH_FASTA_EXPORT
            name = "{0}_{1}.pkl".format(int(G), int(T))
            hjs = os.path.join(folder, name)
            fw = open(hjs, "wb")
            hjs_list.append((folder, name))
        return hjs_list

    def DPI(self):
        for p in os.listdir(self.dp.myCFG.A4_PATH_FASTA_EXPORT):
            if p.endswith(".pkl"):
                os.remove(os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, p))

    def opWdpi(self, list_pep_idx:list):
        fw_list = []
        for pep_idx_info in list_pep_idx:
            hjs = os.path.join(pep_idx_info[0], pep_idx_info[1])
            if os.path.exists(hjs):
                fw = open(hjs, "wb")
                fw_list.append(fw)
            else:
                CLogging.LogForLackPath(hjs)
        return fw_list

    def CloseWriteDistPepIndex(self, fw_list):
        for fw in fw_list:
            fw.close()

    def AddPepIndex(self, pep_info: CPeptideInfo, fw):

        pickle.dump(pep_info, fw)

    def GetPepIndex(self):

        list_pep_index = []
        for p in os.listdir(self.dp.myCFG.A4_PATH_FASTA_EXPORT):
            if not p.endswith('.pkl'): continue
            if p == NAME_PROTEIN_PKL: continue
            if p == NAME_MOD_PKL: continue
            if p.endswith("_ind.pkl"):
                pass
            else:
                assert (os.path.exists(os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, p)))
                assert (os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, p[:-4] + "_ind.pkl"))
                G, T = tool_get_start_T(p)
                list_pep_index.append((p, p[:-4] + "_ind.pkl", G, T))
        list_pep_index.sort(key=lambda x:x[2], reverse=False)
        return list_pep_index