from MSLogging import CLogging

try:
    import cPickle as pickle
except:
    import pickle

def tool_get_start_T(path):
    if path.find('/') >= 0:
        filename = path[path.rfind('/')+1:]
    else:
        filename = path[path.rfind('\\')+1:]
    filename = filename[:filename.find('.')]
    index = filename.find('_')
    G = float(filename[:index])
    T = float(filename[index+1:])
    return G, T

def toolGetWord(inputString, index, d):
    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []

    i = 0
    for c in inputString:
        if c == d:
            p_d.append(i)
        i = i + 1
    result = inputString[p_d[index] + 1:p_d[index + 1]]

    return result

def toolCountCharInString(inputStr, inputChar):

    result = 0

    for c in inputStr:
        if c == inputChar:
            result = result + 1

    return result

def toolGetIndexByWord(line, word, delimiter):

    index = line.find(word)

    if -1 == index:

        CLogging.LogError('can not find: '+word+' from: '+line)

    else:

        return toolCountCharInString(line[0:index], delimiter)

def toolGetNum(inputstr):
    '''
    获取字符串内括号的数字
    :param inputstr: ,str
    :return:数字list
    '''
    index = 0
    left = 0
    right = 0
    list_num = []
    num = 0
    for i in inputstr:
        if (i == '(') | (i == '（'):
            left = index + 1
        elif (i == ')') | (i == '）'):
            right = index
            num = int(inputstr[left: right])
            list_num.append(num)
        index = index + 1
    return list_num

def toolGetListIndex(inputlist, value):
    list_tmp = inputlist.copy()
    out_index = []
    num = list_tmp.count(value)
    for i in range(num):
        replace_index = list_tmp.index(value)
        out_index.append(replace_index)
        list_tmp[replace_index] = -1
    return out_index

def toolCheckLinkSiteValid(C_str, N_str, sq, site, pro_N = False, pro_C = False):
    if not pro_C:
        if site == len(sq)-1 and sq[site] in C_str:
            return False
    if not pro_N:
        if site == 0 and sq[site] in N_str:
            return False
    return True
