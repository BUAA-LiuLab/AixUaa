from MSFlow import *
from MSTask import CTaskG
from MSSystem import EXPIRATION_TIME, INFO_TO_USER_Staff, NAME_FOR_MSREGISTER, SOFTWARE_NAME, LINK_APPLICATION, NAME_PUBLIC_SECURITY
from MSDataC import CDataPack
from MSLogging import CLogging

import datetime
import time
from ctypes import *
import os


class CStaff:
    def __init__(self):
        self.dp = CDataPack()

    def start(self, argv):

        CLogging.LogForStart(SOFTWARE_NAME)
        start_time = time.time()

        self.__captainCheckTime()
        self.__captainRunFlow(argv)

        end_time = time.time()
        CLogging.LogToUser(" Task finish in %.4lfs." % (end_time - start_time))
        CLogging.LogForEnd(SOFTWARE_NAME)

    def __captainCheckTime(self):
        dateNow = datetime.datetime.now()
        dateDead = datetime.datetime(EXPIRATION_TIME['Year'], EXPIRATION_TIME['Month'], EXPIRATION_TIME['Day'], 23, 59)
        n_day = (dateDead - dateNow).days
        if n_day < 0:
            print(INFO_TO_USER_Staff[0])
            exit(-1)
        elif n_day < 7:
            print(INFO_TO_USER_Staff[1])
        else:
            print(INFO_TO_USER_Staff[2])
            print(dateNow)

    def __captainRunFlow(self, argv):
        start_time = time.time()

        if len(argv) == 1:
            flow0 = CFlow0(self.dp)
            flow0.run()

        elif len(argv) == 2:
            print("Config file", argv)
            taskIOConfig = CTaskG(self.dp)
            taskIOConfig.work(argv)
            if self.dp.myCFG.D13_TYPE_TASK == 1:
                flow2 = CFlow2(self.dp)
                flow2.run()
        # 结束
        end_time = time.time()
        print("[Info] Task finish in {}s".format(end_time - start_time))

    def __captainImportDLL(self, path_MSRegister_dll:str=None):

        current_work_dir = os.getcwd()
        if path_MSRegister_dll is None:
            path_MSRegister_dll = os.path.join(current_work_dir, NAME_FOR_MSREGISTER[0])
        if not os.path.exists(path_MSRegister_dll):
            CLogging.LogError("[Error] Missing file : " + NAME_FOR_MSREGISTER[0])
        self.dll = CDLL(path_MSRegister_dll)
        self.dll.RegisterSoftware.argtype = (POINTER(c_char_p))
        self.dll.SetWorkRegister.argtype = (POINTER(c_char_p))
        self.dll.ApplicationLink.argtype = (POINTER(c_char_p))
        self.dll.PublicSecurity.argtype = (POINTER(c_char_p))
        self.dll.RegisterFunction.argtype = (POINTER(c_char_p), POINTER(c_char_p))
        self.dll.RegisterForCheck.restype = (c_bool)
        self.dll.EndOperation.argtype = (c_bool)

    def __captainCheckForRegister(self, path_license:str=None, path_probation:str=None):
        current_work_dir = os.getcwd()
        if path_license is None:
            path_license = os.path.join(current_work_dir, NAME_FOR_MSREGISTER[1])
        if path_probation is None:
            path_probation = os.path.join(current_work_dir, NAME_FOR_MSREGISTER[2])

        software_str = c_char_p(SOFTWARE_NAME.encode())

        application_link = c_char_p(LINK_APPLICATION.encode())
        path_public_security = c_char_p(os.path.join(os.getcwd(), NAME_PUBLIC_SECURITY).encode())

        func1_name = c_char_p("License".encode())
        path1 = c_char_p(path_license.encode())
        func2_name = c_char_p("Probation".encode())
        path2 = c_char_p(path_probation.encode())

        self.dll.RegisterSoftware(software_str)
        self.dll.SetWorkRegister(software_str)
        self.dll.ApplicationLink(application_link)
        self.dll.PublicSecurity(path_public_security)
        self.dll.RegisterFunction(func1_name, path1)
        self.dll.RegisterFunction(func2_name, path2)
        ret = self.dll.RegisterForCheck()
        return ret

    def __captainEndOperation(self, is_exception:bool):

        self.dll.EndOperation(c_bool(is_exception))