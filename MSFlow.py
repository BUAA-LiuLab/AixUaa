from MSTask import CTaskA, CTaskB, CTaskH, CTaskI, CTaskF, CTaskD, CTaskE
from MSTaskSingle import CTaskC, CTaskJ
from MSDataC import CDataPack

from MSSystem import CONFIG_NAME

import os

class CFlow0:
    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP
    def run(self):
        taskA = CTaskA(self.dp)
        taskA.work(os.path.join(os.getcwd(), CONFIG_NAME))


class CFlow1:

    def __init__(self, inputDP:CDataPack):
        self.dp = inputDP

    def run(self):
        taskB = CTaskB(self.dp)
        taskB.work()

        # taskH = CTaskH(self.dp)
        # taskH.work(True)
        if self.dp.myCFG.D16_OPEN_SEARCH_SINGLE != 0:
            taskC = CTaskC(self.dp)
            taskC.search()
            taskJ = CTaskJ(self.dp)
            taskJ.compute_fdr()
        taskI = CTaskI(self.dp)
        taskI.search()
        taskE = CTaskE(self.dp)
        taskE.RerankResult()
        taskD = CTaskD(self.dp)
        taskD.compute_fdr()

# Specific for aa-dependent Uaa cross linking
class CFlow2:
    def __init__(self, inputDP):
        self.dp = inputDP

    def run(self):
        taskB = CTaskB(self.dp)
        taskB.work()

        taskh = CTaskH(self.dp)
        taskh.work(True)
        if self.dp.myCFG.D16_OPEN_SEARCH_SINGLE != 0:
            taskC = CTaskC(self.dp)
            taskC.search()
            taskJ = CTaskJ(self.dp)
            taskJ.compute_fdr()
        taskF = CTaskF(self.dp)
        taskF.search()
        taskE = CTaskE(self.dp)
        taskE.RerankResult()
        taskD = CTaskD(self.dp)
        taskD.compute_fdr()

