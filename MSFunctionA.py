from MSLogging import CLogging
from MSDataC import CConfig

import os
import platform

class CFunctionConfigIO:            # 生成config文件

    @staticmethod
    def config2file(path, config:CConfig):
        with open(path, 'w') as f:
            f.write("#[data]\n")
            f.write("TYPE_MS2=" + config.A1_TYPE_MS2 + "\n")
            f.write("PATH_MS2=" + config.A2_PATH_MS2 + "\n")
            f.write("PATH_FASTA=" + config.A3_PATH_FASTA + "\n")
            f.write("PATH_FASTA_EXPORT=" + config.A4_PATH_FASTA_EXPORT + "\n")
            f.write("PATH_RESULT_EXPORT=" + config.A5_PATH_RESULT_EXPORT + "\n")
            f.write("\n")
            f.write("#[biology]\n")
            f.write("NAME_ENZYME=" + config.B1_NAME_ENZYME + "\n")
            f.write("TYPE_DIGEST=" + config.B2_TYPE_DIGEST + "\n")
            f.write("NUMBER_MAX_MISS_CLV=" + config.B3_NUMBER_MAX_MISS_CLV + "\n")
            f.write("\n")
            f.write("NAME_MOD_FIX=" + config.B4_NAME_MOD_FIX + "\n")
            f.write("NAME_MOD_VAR=" + config.B5_NAME_MOD_VAR + "\n")
            f.write("NUMBER_MAX_MOD=" + config.B6_NUMBER_MAX_MOD + "\n")
            f.write("\n")
            f.write("UAA_SEQ=" + config.B7_UAA_SEQ + "\n")
            f.write("UAA_AA=" + config.B8_UAA_AA + "\n")
            f.write("UAA_LEN_LOW=" + config.B9_UAA_LEN_LOW + "\n")
            f.write("UAA_LEN_UP=" + config.B10_UAA_LEN_UP + "\n")
            f.write("UAA_MASS_LOW=" + config.B11_UAA_MASS_LOW + "\n")
            f.write("UAA_MASS_UP=" + config.B12_UAA_MASS_UP + "\n")
            f.write("UAA_NAME_MOD_FIX=" + config.B13_UAA_NAME_MOD_FIX + "\n")
            f.write("UAA_NAME_MOD_VAR=" + config.B14_UAA_NAME_MOD_VAR + "\n")
            f.write("UAA_COM=" + config.B15_UAA_COM + "\n")
            f.write("UAA_NAME_ENZYME=" + config.B16_UAA_NAME_ENZYME + "\n")
            f.write("UAA_TYPE_DIGEST=" + config.B17_UAA_TYPE_DIGEST + "\n")
            f.write("UAA_NUMBER_MAX_MISS_CLV=" + config.B18_UAA_NUMBER_MAX_MISS_CLV + "\n")
            f.write("UAA_LINKED_AA=" + config.B19_UAA_LINKED_AA + "\n")
            f.write("\n")
            f.write("#[mass spectrometry]\n")
            f.write("TYPE_TOL_PRECURSOR=" + config.C1_TYPE_TOL_PRECURSOR + "\n")
            f.write("PPM_TOL_PRECURSOR=" + config.C2_PPM_TOL_PRECURSOR + "\n")
            f.write("TYPE_TOL_FRAGMENT=" + config.C3_TYPE_TOL_FRAGMENT + "\n")
            f.write("PPM_TOL_FRAGMENT=" + config.C4_PPM_TOL_FRAGMENT + "\n")
            f.write("TYPE_ACTIVATION=" + config.C5_ACTIVATION_TYPE + "\n")
            f.write("\n")
            f.write("#[performance]\n")
            f.write("NUMBER_THREAD=" + config.D1_NUMBER_THREAD + "\n")
            f.write("TYPE_THREAD=" + config.D2_TYPE_THREAD + "\n")
            f.write("NUMBER_SELECT_PEAK=" + config.D3_NUMBER_SELECT_PEAK + "\n")
            f.write("NUMBER_SPECTRUM=" + config.D4_NUMBER_SPECTRUM + "\n")
            f.write("LEN_MAX_PROTEIN=" + config.D5_LEN_MAX_PROTEIN + "\n")
            f.write("MASS_PEP_LOW=" + config.D6_MASS_PEP_LOW + "\n")
            f.write("MASS_PEP_UP=" + config.D7_MASS_PEP_UP + "\n")
            f.write("LEN_PEP_LOW=" + config.D8_LEN_PEP_LOW + "\n")
            f.write("LEN_PEP_UP=" + config.D9_LEN_PEP_UP + "\n")
            f.write("INDEX_SPLIT_MASS=" + config.D10_INDEX_SPLIT_MASS + "\n")
            f.write("NUMBER_TOP_RESULT=" + config.D11_NUMBER_TOP_RESULT + "\n")
            f.write("\n")
            f.write("MULTI_MASS=" + config.D12_MULTI_MASS + "\n")
            f.write("TYPE_TASK=" + config.D13_TYPE_TASK + "\n")
            f.write("TYPE_FILTER_BETA=" + config.D14_TYPE_FILTER_BETA + "\n")
            f.write("NUMBER_PEAK_BETA=" + config.D15_NUMBER_PEAK_BETA + "\n")
            f.write("\n")
            f.write("OPEN_SEARCH_SINGLE=" + config.D16_OPEN_SEARCH_SINGLE + "\n")
            f.write("MASS_WINDOW_BETA=" + config.D17_MASS_WINDOW_BETA + "\n")
            f.write("PATH_PFIND_RESULT=" + config.D18_PATH_PFIND_RESULT + "\n")
            f.write("\n")
            f.write("#[filter]\n")
            f.write("FDR_PSM=" + config.E1_FDR_PSM + "\n")
            f.write("\n")
            f.write("#[ini]\n")
            f.write("PATH_INI_ELEMENT=" + config.F1_PATH_INI_ELEMENT + "\n")
            f.write("PATH_INI_AA=" + config.F2_PATH_INI_AA + "\n")
            f.write("PATH_INI_MOD=" + config.F3_PATH_INI_MOD + "\n")
            f.write("\n")
            f.write("#[advance]\n")
            f.write("CLV_UAA_COM=" + config.G2_BETA_UAA_NEW_COM + "\n")
            f.write("CLV_UAA_LINKED_AA=" + config.G3_ALPHA_AA_TYPE + "\n")
            f.write("CLV_UAA_LINKED_AA_MASS_CHANGE=" + config.G4_ALPHA_AA_MASS_CHANGE + "\n")

    @staticmethod
    def file2config(file_path, config:CConfig):
        sys_flag = platform.system()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.find('=') < 0:
                    continue
                if line.find('#') >= 0:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2:
                    continue
                name = tmp[0]
                value = tmp[1]
                print(name, ":", value)
                try:
                    # [ini]
                    if name == "PATH_INI_ELEMENT":
                        config.F1_PATH_INI_ELEMENT = value
                    elif name == "PATH_INI_AA":
                        try:
                            config.F2_PATH_INI_AA = value
                        except Exception as err:
                            print("[Error] The element_file must be in front of AA_file.", err)
                            info = "[Error] The element_file must be in front of AA_file."
                            CLogging.LogError(info)
                    elif name == "PATH_INI_MOD":
                        config.F3_PATH_INI_MOD = value

                    # [data]
                    elif name == "TYPE_MS2":
                        config.A1_TYPE_MS2 = value
                    elif name == "PATH_MS2":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A2_PATH_MS2 = value
                    elif name == "PATH_FASTA":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A3_PATH_FASTA = value
                    elif name == "PATH_FASTA_EXPORT":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A4_PATH_FASTA_EXPORT = value
                        if not os.path.exists(config.A4_PATH_FASTA_EXPORT):
                            os.makedirs(config.A4_PATH_FASTA_EXPORT)
                    elif name == "PATH_RESULT_EXPORT":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A5_PATH_RESULT_EXPORT = value
                        if not os.path.exists(config.A5_PATH_RESULT_EXPORT):
                            os.makedirs(config.A5_PATH_RESULT_EXPORT)

                    # [biology]
                    elif name == "NAME_ENZYME":
                        config.B1_NAME_ENZYME = value
                    elif name == "TYPE_DIGEST":
                        config.B2_TYPE_DIGEST = int(value)
                    elif name == "NUMBER_MAX_MISS_CLV":
                        config.B3_NUMBER_MAX_MISS_CLV = int(value)
                    elif name == "NAME_MOD_FIX":
                        config.B4_NAME_MOD_FIX = value
                    elif name == "NAME_MOD_VAR":
                        config.B5_NAME_MOD_VAR = value
                    elif name == "NUMBER_MAX_MOD":
                        config.B6_NUMBER_MAX_MOD = int(value)
                    elif name == "UAA_SEQ":
                        config.B7_UAA_SEQ = value
                    elif name == "UAA_AA":
                        config.B8_UAA_AA  = value[0]
                    elif name == "UAA_LEN_LOW":
                        config.B9_UAA_LEN_LOW = int(value)
                    elif name == "UAA_LEN_UP":
                        config.B10_UAA_LEN_UP = int(value)
                    elif name == "UAA_MASS_LOW":
                        config.B11_UAA_MASS_LOW = float(value)
                    elif name == "UAA_MASS_UP":
                        config.B12_UAA_MASS_UP = float(value)
                    elif name == "UAA_NAME_MOD_FIX":
                        config.B13_UAA_NAME_MOD_FIX = value
                    elif name == "UAA_NAME_MOD_VAR":
                        config.B14_UAA_NAME_MOD_VAR = value
                    elif name == "UAA_COM":
                        config.B15_UAA_COM = value
                    elif name == "UAA_NAME_ENZYME":
                        config.B16_UAA_NAME_ENZYME = value
                    elif name == "UAA_TYPE_DIGEST":
                        config.B17_UAA_TYPE_DIGEST = int(value)
                    elif name == "UAA_NUMBER_MAX_MISS_CLV":
                        config.B18_UAA_NUMBER_MAX_MISS_CLV = int(value)
                    elif name == "UAA_LINKED_AA":
                        config.B19_UAA_LINKED_AA = value

                    # [mass spectrometry]
                    elif name == "TYPE_TOL_PRECURSOR":
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.C1_PPM_TOL_PRECURSOR = True
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.C1_PPM_TOL_PRECURSOR = False
                        else:
                            CLogging.LogError(line)

                    elif name == "PPM_TOL_PRECURSOR":
                        if config.C1_PPM_TOL_PRECURSOR:
                            config.C2_PPM_TOL_PRECURSOR = float(value) * 1e-6
                        else:
                            config.C2_PPM_TOL_PRECURSOR = float(value)  # default is 20ppm

                    elif name == "TYPE_TOL_FRAGMENT":
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.C3_TYPE_TOL_FRAGMENT = True
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.C3_TYPE_TOL_FRAGMENT = False
                        else:
                            CLogging.LogError(line)

                    elif name == "PPM_TOL_FRAGMENT":
                        if config.C3_TYPE_TOL_FRAGMENT:
                            config.C4_PPM_TOL_FRAGMENT = float(value) * 1e-6
                        else:
                            config.C4_PPM_TOL_FRAGMENT = float(value)  # default is 20ppm

                    elif name == "PPM_TOL_FRAGMENT":
                        config.is_ppm_fra = True
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.C4_PPM_TOL_FRAGMENT = float(value[:-3]) * 1e-6
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.is_ppm_fra = False
                            config.C4_PPM_TOL_FRAGMENT = float(value[:-2])
                        else:
                            config.C4_PPM_TOL_FRAGMENT = 20e-6  # default is 20ppm

                    # [performance]
                    elif name == "NUMBER_THREAD":
                        config.D1_NUMBER_THREAD = int(value)
                    elif name == "TYPE_THREAD":
                        config.D2_TYPE_THREAD = int(value)
                    elif name == "NUMBER_SELECT_PEAK":
                        config.D3_NUMBER_SELECT_PEAK = int(value)
                    elif name == "NUMBER_SPECTRUM":
                        config.D4_NUMBER_SPECTRUM = int(value)

                    elif name == "MASS_PEP_LOW":
                        config.D6_MASS_PEP_LOW  = float(value)
                    elif name == "MASS_PEP_UP":
                        config.D7_MASS_PEP_UP = float(value)
                    elif name == "LEN_PEP_LOW":
                        config.D8_LEN_PEP_LOW = int(value)
                    elif name == "LEN_PEP_UP":
                        config.D9_LEN_PEP_UP = int(value)
                    elif name == "INDEX_SPLIT_MASS":
                        config.D10_INDEX_SPLIT_MASS = float(value)
                    elif name == "NUMBER_TOP_RESULT":
                        config.D11_NUMBER_TOP_RESULT = int(value)
                    elif name == "MULTI_MASS":
                        config.D12_MULTI_MASS = int(value)
                    elif name == "TYPE_TASK":
                        config.D13_TYPE_TASK = int(value)
                    elif name == "TYPE_FILTER_BETA":
                        config.D14_TYPE_FILTER_BETA = (int(value) == 1)
                    elif name == "NUMBER_PEAK_BETA":
                        config.D15_NUMBER_PEAK_BETA = int(value)
                    elif name == "OPEN_SEARCH_SINGLE":
                        config.D16_OPEN_SEARCH_SINGLE = int(value)
                    elif name == "MASS_WINDOW_BETA":
                        config.D17_MASS_WINDOW_BETA  = float(value)
                    elif name == "PATH_PFIND_RESULT":
                        config.D18_PATH_PFIND_RESULT = value.strip()

                    # [filter]
                    elif name == "FDR_PSM":
                        config.E1_FDR_PSM = float(value)
                        config.using_decoy = True

                    # [advance]
                    elif name == "CLV_UAA_COM":
                        config.G2_BETA_UAA_NEW_COM = value
                    elif name == "CLV_UAA_LINKED_AA":
                        config.G3_ALPHA_AA_TYPE = value
                    elif name == "CLV_UAA_LINKED_AA_MASS_CHANGE":
                        config.G4_ALPHA_AA_MASS_CHANGE = float(value)

                except Exception as err:
                    CLogging.LogError(f"[Error] occured in {file_path}" + str(err))
