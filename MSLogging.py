# -*- mode: python ; coding: utf-8 -*-
import logging
import sys
import traceback
import time

from MSSystem import SOFTWARE_NAME

PATH_LOG = f'{SOFTWARE_NAME}.log'
PATH_ERROR = f'{SOFTWARE_NAME}_error.log'


def Assert(value:bool):
    if not value:
        tb_stack = traceback.extract_stack()
        info = "[Error] "
        for s in tb_stack:
            info += s.filename + ' ' + str(s.lineno) + ' ' + s.name + '\n'
        CLogging.LogError(info)
        assert value

class CLogging:

    @staticmethod
    def GetLocTimeStr():
        loc_time = time.localtime()
        str_loc_time = str(loc_time.tm_year) + "-" + str(loc_time.tm_mon) + "-" + str(loc_time.tm_mday) + "-" + str(
            loc_time.tm_hour) + "-" + str(loc_time.tm_min) + "-" + str(loc_time.tm_sec)
        return str_loc_time

    @staticmethod
    def LogForStart(str_name:str=""):
        strInfo = "=========== Start " + str_name + " ==========="
        CLogging.LogToUser(strInfo)

    @staticmethod
    def LogForEnd(str_name:str=""):
        strInfo = "=========== End " + str_name + " ==========="
        CLogging.LogToUser(strInfo)

    @staticmethod
    def LogForLackPath(path:str, info_type:int=1):

        warn_str_name = "Path is not exist ! " + path
        if info_type == 1:
            CLogging.LogError(warn_str_name)
        elif info_type == 2:
            CLogging.LogWarning(warn_str_name)
        else:
            CLogging.LogError(warn_str_name)

    @staticmethod
    def LogToUser(strInfo:str, log2file:bool=True, log2user:bool=True, log_type:str="Info"):
        Info = "[" + log_type + "] " + CLogging.GetLocTimeStr() + " " + strInfo
        if log2user:
            print(Info)
        # if os.access(myLogPath, os.Ws_OK):  # 当文件被excel打开时，这个东东没啥用
        if log2file:
            try:
                f_w = open(PATH_LOG, 'a', encoding='utf8')
                f_w.write(Info + '\n')
                f_w.close()
            except IOError:
                print("CodeTemplate.log is opened! Please close it and run the program again!")
                sys.exit(0)

    @staticmethod
    def LogError(info):


        logging.basicConfig(filename=PATH_LOG,
                            filemode='a',
                            format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                            )

        logging.basicConfig(filename=PATH_ERROR,
                            filemode='a',
                            format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                            )

        stack_info = str(info) + "\n"
        exc_type, exc_value, exc_trace_back = sys.exc_info()
        while exc_trace_back is not None:
            stack_info += exc_trace_back.tb_frame.f_globals["__file__"] + ' ' + str(exc_trace_back.tb_lineno) + '\n'
            exc_trace_back = exc_trace_back.tb_next
        CLogging.LogToUser(stack_info, log_type="Error")
        logging.error(stack_info)
        sys.exit("Error")
        # sys.exit(0)

    @staticmethod
    def LogWarning(info):
        CLogging.LogToUser(info, log_type="Warn")

        logging.basicConfig(filename=PATH_LOG,
                            filemode='a',
                            format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                            )
        logging.warning(info)