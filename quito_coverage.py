from qiskit import (
    #IBMQ,
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
    execute,
    Aer,
)

import numpy as np
import os
import sys

import configparser
import os.path
import importlib
import time
import warnings
from scipy.stats import wilcoxon
# from qiskit.providers.aer import AerError

def quito(
        con_file: str
):
    warnings.filterwarnings('ignore')
    _quito_run(con_file)

def _end_running():
    exit()

def _check_configuration_file(config):
    if config.has_section('program') == False:
        print("Error: Quito cannot find section 'program' in this configuration file.")
        _end_running()
    else:
        if config.has_option('program', 'root') == False:
            print("Error: Quito cannot find the root of the program.")
            _end_running()
        if config.has_option('program', 'num_qubit') == False:
            print("Error: Quito cannot find the number of qubits of the program.")
            _end_running()
        if config.has_option('program', 'inputID') == False:
            print("Error: Quito cannot find the input IDs of the program.")
            _end_running()
        if config.has_option('program', 'outputID') == False:
            print("Error: Quito cannot find the output IDs of the program.")
            _end_running()
    if config.has_section('program_specification_category') == False:
        print("Error: Quito cannot find section 'program_specification_category' in this configuration file.")
        _end_running()
    else:
        if config.has_option('program_specification_category', 'ps_category') == False:
            print("Error: Quito cannot find the category of the program specification.")
            _end_running()
    if config.has_section('quito_configuration') == False:
        print("Error: Quito cannot find section 'quito_configuration' in this configuration file.")
        _end_running()
    else:
        if config.has_option('quito_configuration', 'coverage_criterion') == False:
            print("Error: Quito cannot find the coverage criterion you choose.")
            _end_running()

    ps_category = config.get('program_specification_category', 'ps_category')
    if ps_category == 'full' or ps_category == 'partial':
        if config.has_section('program_specification') == False:
            print("Error: Quito cannot find the program specification.")
            _end_running()
    return ps_category

def _check_inputID_outputID(num_qubit, inputID, outputID):
    if _check_unique(inputID) == False:
        print("Wrong input IDs")
        _end_running()
    if _check_unique(outputID) == False:
        print("Wrong output IDs")
        _end_running()
    inputID.sort()
    outputID.sort()

    if int(inputID[-1]) > num_qubit - 1:
        print("Wrong input IDs")
        _end_running()
    if int(inputID[-1]) > num_qubit - 1:
        print("Wrong output IDs")
        _end_running()

    return inputID, outputID

def _check_unique(l):
    return len(l) == len(set(l))


def _check_same(l,value):
    for i in range(len(l)):
        if l[i] != value:
            return False
    return True


def _check_WOO(i, o , valid_inputs, valid_outputs):
    flag = False
    for k in range(len(valid_inputs)):
        if valid_inputs[k] == i and valid_outputs[k] == o:
            flag = True
    if flag == False:
        print('fail for woo')
        return False
    return True


def _wilcoxon(fre,p):
    pvalue = []
    res = []
    for i in range(len(p)):
        if np.isnan(fre[i]).any() == True:
            pvalue.append(-1)
            res.append(-1)
        elif _check_same(fre[i],p[i]) == True:
            pvalue.append(1)
            res.append(1)
        else:
            fre_np = np.array(fre[i], dtype=float)
            result = wilcoxon(fre_np, correction=True, y=np.repeat(p[i], len(fre[i])))
            pvalue.append(result[1])
    return pvalue

def _judge_ass_result(inputs, outputs, pvalue, C_LEVEL, f):
    for i in range(len(pvalue)):
        if pvalue[i] == -1:
            f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + "Not enough inputs for statistical test.")
            f.write('\n')
        elif pvalue[i] < C_LEVEL:
            f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + str(pvalue[i]) + "; ")
            f.write("Result: Fail for OPO")
            f.write('\n')
        elif pvalue[i] >= C_LEVEL:
            f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + str(pvalue[i]) + "; ")
            f.write("Result: Inconclusive")
            f.write('\n')

# def _process_bar(percent, start_str='', end_str='', total_length=0):
#     bar = ''.join(['#'] * int(percent * total_length)) + ''
#     bar = '\r' + start_str + bar.ljust(total_length) + ' {:0>4.1f}%|'.format(percent*100) + end_str
#     print(bar, end='', flush=True)

def _check_partial_ps(valid_input, valid_output, p):#check whether there is any input having full ps(sum of p of all outputs of one input is 1)
    index = _input_group(valid_input)
    p_sum = 0
    input_full = []
    output_full = []
    for i in range(len(index)):
        start = index[i]
        if i == len(index) - 1:
            end = len(valid_input)
        else:
            end = index[i+1]
        for j in range(start,end):
            p_sum += p[j]
        if p_sum == 1:
            for k in range(start,end):
                input_full.append(valid_input[k])
                output_full.append(valid_output[k])
        p_sum = 0
    return input_full, output_full

def _get_unique(l):
    unique = []
    for i in l:
        if i not in unique:
            unique.append(i)
    return unique

def _dec2bin(n,bit):
    a = 1
    list = []
    while a > 0:
        a, b = divmod(n, 2)
        list.append(str(b))
        n = a
    s = ""
    for i in range(len(list) - 1, -1, -1):
        s += str(list[i])
    s = s.zfill(bit)
    return s

def _get_all(bit):
    all = []
    for i in range(pow(2,bit)):
        i_bin = _dec2bin(i, bit)
        all.append(i_bin)
    return all

def _ps_dic(valid_input, valid_output):
    dic = {}
    input_flag = valid_input[0]
    dic[input_flag] = [valid_output[0]]
    for i in range(1, len(valid_input)):
        if valid_input[i] == input_flag:
            dic[input_flag].append(valid_output[i])
        else:
            input_flag = valid_input[i]
            dic[input_flag] = [valid_output[i]]
    return dic

def _com_dic(all_input, all_output):
    dic = {}
    for input in all_input:
        dic[input] = []
        for output in all_output:
            dic[input].append(output)
    return dic

def _execute_quantum_program(inputID, outputID, num_qubit, i, module_name, shots):
    q = QuantumRegister(num_qubit)
    c = ClassicalRegister(len(outputID))
    qc = QuantumCircuit(q, c)
    for j in range(len(inputID)):
        if i[len(inputID) - 1 - j] == '1':
            qc.x(int(inputID[j]))
    module = importlib.import_module(module_name)
    run_method = getattr(module,"run")
    run_method(qc)
    simulator = Aer.get_backend('aer_simulator')
    simulator.set_options(device='GPU')
    result = execute(qc, simulator, shots=shots, memory=True).result().get_memory()

    return result

def _input_coverage(con_dic):
    unique_inputs = _get_unique(con_dic['valid_input']) #unique inputs list
    r = int(con_dic['K'] / con_dic['M'])  # group size
    counts = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    result_dic = {}
    input_data = []
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    T2 = 0
    for i_index in range(len(unique_inputs)):
        i = unique_inputs[i_index]
        T2_start = time.time()
        # Each output value in memory list is for one test suite
        memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                          con_dic['module_name'], con_dic['K'])
        T2_end = time.time()
        T2 += T2_end - T2_start
        input_data.append(unique_inputs[i_index])
        #Check WOO
        for o_index in range(len(memory)):
            o = memory[o_index]
            group_num = int(o_index/r)
            result_dic['output_data'][o_index].append(o)
            if not _check_WOO(i, o, con_dic['valid_input'], con_dic['valid_output']): #fail for woo
                result_dic['WOO'] = True
                result_dic['w_input'] = {i:o}
                result_dic['input_data'] = [input_data]*con_dic['K']
                result_dic['T2'] = T2
                return result_dic
            for mark in range(len(con_dic['valid_input'])):
                if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                    counts[mark][group_num] += 1
    # If there is no unexpected outputs
    fre = [[e/r for e in row] for row in counts]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['input_data'] = [input_data]*con_dic['K']
    result_dic['T2'] = T2
    return result_dic

def _input_coverage_partial(con_dic):
    # getting all possible inputs
    all_inputs = _get_all(len(con_dic['inputID']))
    r = int(con_dic['K'] / con_dic['M'])  # group size
    counts = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    result_dic = {}
    input_data = []
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    T2 = 0

    # Check whether there is full ps for some input, for those we need to check WOO
    input_full, output_full = _check_partial_ps(con_dic['valid_input'], con_dic['valid_output'], con_dic['p'])

    for i_index in range(len(all_inputs)):
        i = all_inputs[i_index]
        T2_start = time.time()
        memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                          con_dic['module_name'], con_dic['K'])
        T2_end = time.time()
        T2 += T2_end - T2_start
        input_data.append(i)
        for o_index in range(len(memory)):
            o = memory[o_index]
            group_num = int(o_index/r)
            result_dic['output_data'][o_index].append(o)
            if i in input_full:
                # Check WOO
                if not _check_WOO(i,o,input_full,output_full):
                    result_dic['WOO'] = True
                    result_dic['input_data'] = [input_data] * con_dic['K']
                    result_dic['w_input'] = {i: o}
                    result_dic['T2'] = T2
                    return result_dic
            for mark in range(len(con_dic['valid_input'])):
                if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                    counts[mark][group_num] += 1
    # If there is no unexpected outputs
    fre = [[e / r for e in row] for row in counts]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['input_data'] = [input_data] * con_dic['K']
    result_dic['T2'] = T2
    return result_dic



#
# def input_coverage_partial(con_dic):
#     resultFolder = os.path.join(con_dic['program_folder'], "result")
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = _get_all(len(con_dic['inputID']))
#     #input_num = len(unique_inputs)
#     r = int(con_dic['K']/con_dic['M'])
#     counts = np.zeros((len(con_dic['valid_input']),con_dic['M']))
#     fre = np.zeros((len(con_dic['valid_input']),con_dic['M']))
#     co = 0
#     count_cases = 0
#     sum = np.zeros((len(all_inputs),con_dic['M']))
#     input_full, output_full = _check_partial_ps(con_dic['valid_input'], con_dic['valid_output'], con_dic['p'])
#     # print("input full:"+str(input_full))
#     # print("output full:"+str(output_full))
#
#     input_file = open(resultFolder + 'INPUTS_input_coverage_'+con_dic['module_name']+'.txt','w')
#     result_file = open(resultFolder + 'RESULTS_input_coverage_'+con_dic['module_name']+'.txt','w')
#     ass_file = open(resultFolder + 'ASSESSMENT_input_coverage_'+con_dic['module_name']+'.txt','w')
#     T2 = 0
#     T1_start = time.time()
#     for g in range(con_dic['M']):
#         for t in range(r):
#             co += 1
#             _process_bar(co / con_dic['K'], start_str='', end_str='100%', total_length=30)
#             for i in all_inputs: #i为input
#                 count_cases += 1
#                 start = time.time()
#                 result = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i, con_dic['module_name'])
#                 end = time.time()
#                 T2 += end - start
#                 input_file.write(str(i) + ' ')
#                 sum[int(i,2)][g] += 1
#                 # print(i)
#                 # print(result)
#                 o = list(result)[0]
#                 result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                 if i in input_full:
#                     if _check_WOO(i, o, input_full, output_full) == False:
#                         T1_end = time.time()
#                         T1 = T1_end - T1_start
#                         ass_file.write("Fail for WOO")
#                         result_file.write('\n')
#                         result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                         result_file.write('\n')
#                         result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                         print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         _end_running()
#                 #print('input='+i)
#                 #print('output='+o)
#                 for mark in range(len(con_dic['valid_input'])):
#                     if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
#                         counts[mark][g] += 1
#             input_file.write('\n')
#             result_file.write('\n')
#     print('\n')
#
#     T1_end = time.time()
#     T1 = T1_end - T1_start
#
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2)+"s.")
#
#     #getting frequency
#     #sum =np.zeros(M)
#     # print(counts)
#     # print(sum)
#     input_index = _input_group(con_dic['valid_input'])
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(con_dic['valid_input'])
#         else:
#             end = input_index[i+1]
#         for j in range(start, end):
#             fre[j] = counts[j]/sum[int(con_dic['valid_input'][start],2)]
#         # print(sum)
#
#     # print(fre)
#     pvalue = _wilcoxon(fre, con_dic['p'])
#     _judge_ass_result(con_dic['valid_input'], con_dic['valid_output'], pvalue, ass_file)
#
#     print("The total run time is " + "{:.2f}".format(T1) + "s.")
#     print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()

def _input_coverage_no(con_dic):
    # getting all possible inputs
    all_inputs = _get_all(len(con_dic['inputID']))
    result_dic = {}
    input_data = []
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    T2 = 0

    for i_index in range(len(all_inputs)):
        i = all_inputs[i_index]
        T2_start = time.time()
        memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                          con_dic['module_name'], con_dic['K'])
        T2_end = time.time()
        T2 += T2_end - T2_start
        input_data.append(i)
        for o_index in range(len(memory)):
            o = memory[o_index]
            result_dic['output_data'][o_index].append(o)
    result_dic['input_data'] = [input_data] * con_dic['K']
    result_dic['T2'] = T2
    return result_dic

def _output_coverage(con_dic):
    unique_inputs = _get_unique(con_dic['valid_input']) #unique inputs list
    unique_outputs = _get_unique(con_dic['valid_output']) #unique outputs list
    output_flag = np.zeros([con_dic['K'], len(unique_outputs)])
    r = int(con_dic['K'] / con_dic['M'])  # group size
    counts = np.zeros([len(con_dic['valid_input']), con_dic['M']])
    counts_input = np.zeros([len(unique_inputs), con_dic['M']])
    result_dic = {}
    # input_data = []
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    shots = con_dic['K']
    exe_count = 0
    stop_flag = False
    not_finish_list = [i for i in range(con_dic['K'])] #Test suites that haven't got all output values
    T2 = 0
    #The maximum test suite size is the budget
    while True:
        for i_index in range(len(unique_inputs)):
            exe_count += 1
            i = unique_inputs[i_index]
            pop_list = []
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                          con_dic['module_name'], shots)
            T2_end = time.time()
            T2 += T2_end-T2_start
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                group_num = int(not_finish_list[o_index]/r)
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                if not _check_WOO(i, o, con_dic['valid_input'], con_dic['valid_output']):  # fail for woo
                    result_dic['WOO'] = True
                    result_dic['w_input'] = {i: o}
                    result_dic['T2'] = T2
                    return result_dic
                for uo_index in range(len(unique_outputs)):
                    uo = unique_outputs[uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]],1):
                    pop_list.append(o_index)
                    # not_finish_list.pop(o_index)
                    shots -= 1
                for mark in range(len(con_dic['valid_input'])):
                    if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                        counts[mark][group_num] += 1
                        counts_input[unique_inputs.index(i)][group_num] += 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                stop_flag = True
                break
            if exe_count == con_dic['BUDGET']:
                stop_flag = True
                break
        if stop_flag:
            break
    fre = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    for i_index in range(len(con_dic['valid_input'])):
        for g_index in range(con_dic['M']):
            i = con_dic['valid_input'][i_index]
            fre[i_index][g_index] = counts[i_index][g_index]/counts_input[unique_inputs.index(i)][g_index]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['T2'] = T2
    return result_dic

def _output_coverage_partial(con_dic):
    all_inputs = _get_all(len(con_dic['inputID']))
    all_outputs = _get_all(len(con_dic['outputID']))
    valid_unique_inputs = _get_unique(con_dic['valid_input'])
    output_flag = np.zeros([con_dic['K'], len(all_outputs)])
    r = int(con_dic['K'] / con_dic['M'])  # group size
    counts = np.zeros([len(con_dic['valid_input']), con_dic['M']])
    counts_input = np.zeros([len(valid_unique_inputs), con_dic['M']])
    result_dic = {}
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    shots = con_dic['K']
    exe_count = 0
    stop_flag = False
    not_finish_list = [i for i in range(con_dic['K'])]  # Test suites that haven't got all output values

    # Check whether there is full ps for some input, for those we need to check WOO
    input_full, output_full = _check_partial_ps(con_dic['valid_input'], con_dic['valid_output'], con_dic['p'])
    T2 = 0
    # The maximum test suite size is the budget
    while True:
        for i_index in range(len(all_inputs)):
            exe_count += 1
            i = all_inputs[i_index]
            pop_list = []
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                              con_dic['module_name'], shots)
            T2_end = time.time()
            T2 += T2_end - T2_start
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                group_num = int(not_finish_list[o_index] / r)
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                if i in input_full:
                    if not _check_WOO(i, o, input_full, output_full):  # fail for woo
                        result_dic['WOO'] = True
                        result_dic['w_input'] = {i: o}
                        result_dic['T2'] = T2
                        return result_dic
                for uo_index in range(len(all_outputs)):
                    uo = all_outputs[uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]], 1):
                    pop_list.append(o_index)
                    shots -= 1
                for mark in range(len(con_dic['valid_input'])):
                    if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                        counts[mark][group_num] += 1
                        counts_input[valid_unique_inputs.index(i)][group_num] += 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                stop_flag = True
                break
            if exe_count == con_dic['BUDGET']:
                stop_flag = True
                break
        if stop_flag:
            break
    fre = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    for i_index in range(len(con_dic['valid_input'])):
        for g_index in range(con_dic['M']):
            i = con_dic['valid_input'][i_index]
            fre[i_index][g_index] = counts[i_index][g_index] / counts_input[valid_unique_inputs.index(i)][g_index]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['T2'] = T2
    return result_dic

def _output_coverage_no(con_dic):
    all_inputs = _get_all(len(con_dic['inputID']))
    all_outputs = _get_all(len(con_dic['outputID']))
    output_flag = np.zeros([con_dic['K'], len(all_outputs)])
    result_dic = {}
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    shots = con_dic['K']
    exe_count = 0
    stop_flag = False
    not_finish_list = [i for i in range(con_dic['K'])]  # Test suites that haven't got all output values
    T2 = 0

    # The maximum test suite size is the budget
    while True:
        for i_index in range(len(all_inputs)):
            exe_count += 1
            i = all_inputs[i_index]
            pop_list = []
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                              con_dic['module_name'], shots)
            T2_end = time.time()
            T2 += T2_end - T2_start
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                for uo_index in range(len(all_outputs)):
                    uo = all_outputs[uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]], 1):
                    pop_list.append(o_index)
                    shots -= 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                stop_flag = True
                break
            if exe_count == con_dic['BUDGET']:
                stop_flag = True
                break
        if stop_flag:
            break
    result_dic['T2'] = T2
    return result_dic

def _input_output_coverage(con_dic):
    unique_inputs = _get_unique(con_dic['valid_input'])  # unique inputs list
    BUDGET_I = int(con_dic['BUDGET'] / len(unique_inputs)) # budget for each input
    ps_dic = _ps_dic(con_dic['valid_input'], con_dic['valid_output'])
    r = int(con_dic['K'] / con_dic['M'])  # group size
    result_dic = {}
    counts = np.zeros([len(con_dic['valid_input']), con_dic['M']])
    counts_input = np.zeros([len(unique_inputs), con_dic['M']])
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    T2 = 0

    for i in ps_dic.keys():
        exe_count = 0
        not_finish_list = [i for i in range(con_dic['K'])]
        output_flag = np.zeros([con_dic['K'], len(ps_dic[i])])
        shots = con_dic['K']
        while True:
            exe_count += 1
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                                      con_dic['module_name'], shots)
            T2_end = time.time()
            T2 += T2_end - T2_start
            pop_list = []
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                group_num = int(not_finish_list[o_index]/r)
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                if not _check_WOO(i, o, con_dic['valid_input'], con_dic['valid_output']):  # fail for woo
                    result_dic['WOO'] = True
                    result_dic['w_input'] = {i: o}
                    result_dic['T2'] = T2
                    return result_dic
                for uo_index in range(len(ps_dic[i])):
                    uo = ps_dic[i][uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]], 1):
                    pop_list.append(o_index)
                    # not_finish_list.pop(o_index)
                    shots -= 1
                for mark in range(len(con_dic['valid_input'])):
                    if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                        counts[mark][group_num] += 1
                        counts_input[unique_inputs.index(i)][group_num] += 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                break
            if exe_count == BUDGET_I:
                break
    fre = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    for i_index in range(len(con_dic['valid_input'])):
        for g_index in range(con_dic['M']):
            i = con_dic['valid_input'][i_index]
            fre[i_index][g_index] = counts[i_index][g_index] / counts_input[unique_inputs.index(i)][g_index]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['T2'] = T2
    return result_dic

def _input_output_coverage_partial(con_dic):
    all_inputs = _get_all(len(con_dic['inputID']))  # unique inputs list
    all_outputs = _get_all(len(con_dic['outputID']))
    valid_unique_inputs = _get_unique(con_dic['valid_input'])
    BUDGET_I = int(con_dic['BUDGET'] / len(all_inputs)) # budget for each input
    com_dic = _com_dic(all_inputs, all_outputs)
    r = int(con_dic['K'] / con_dic['M'])  # group size
    result_dic = {}
    counts = np.zeros([len(con_dic['valid_input']), con_dic['M']])
    counts_input = np.zeros([len(valid_unique_inputs), con_dic['M']])
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['WOO'] = False
    T2 = 0

    input_full, output_full = _check_partial_ps(con_dic['valid_input'], con_dic['valid_output'], con_dic['p'])

    for i in com_dic.keys():
        exe_count = 0
        not_finish_list = [i for i in range(con_dic['K'])]
        output_flag = np.zeros([con_dic['K'], len(com_dic[i])])
        shots = con_dic['K']
        while True:
            exe_count += 1
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                                      con_dic['module_name'], shots)
            T2_end = time.time()
            T2 += T2_end - T2_start
            pop_list = []
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                group_num = int(not_finish_list[o_index] / r)
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                if i in input_full:
                    if not _check_WOO(i, o, input_full, output_full):  # fail for woo
                        result_dic['WOO'] = True
                        result_dic['w_input'] = {i: o}
                        result_dic['T2'] = T2
                        return result_dic
                for uo_index in range(len(com_dic[i])):
                    uo = com_dic[i][uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]], 1):
                    pop_list.append(o_index)
                    # not_finish_list.pop(o_index)
                    shots -= 1
                for mark in range(len(con_dic['valid_input'])):
                    if con_dic['valid_input'][mark] == i and con_dic['valid_output'][mark] == o:
                        counts[mark][group_num] += 1
                        counts_input[valid_unique_inputs.index(i)][group_num] += 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                break
            if exe_count == BUDGET_I:
                break
    fre = np.zeros((len(con_dic['valid_input']), con_dic['M']))
    for i_index in range(len(con_dic['valid_input'])):
        for g_index in range(con_dic['M']):
            i = con_dic['valid_input'][i_index]
            fre[i_index][g_index] = counts[i_index][g_index] / counts_input[valid_unique_inputs.index(i)][g_index]
    result_dic['pvalue'] = _wilcoxon(fre, con_dic['p'])
    result_dic['T2'] = T2
    return result_dic

def _input_output_coverage_no(con_dic):
    all_inputs = _get_all(len(con_dic['inputID']))  # unique inputs list
    all_outputs = _get_all(len(con_dic['outputID']))
    BUDGET_I = int(con_dic['BUDGET'] / len(all_inputs)) # budget for each input
    com_dic = _com_dic(all_inputs, all_outputs)
    result_dic = {}
    result_dic['input_data'] = [[] for _ in range(con_dic['K'])]
    result_dic['output_data'] = [[] for _ in range(con_dic['K'])]
    T2 = 0

    for i in com_dic.keys():
        exe_count = 0
        not_finish_list = [i for i in range(con_dic['K'])]
        output_flag = np.zeros([con_dic['K'], len(com_dic[i])])
        shots = con_dic['K']
        while True:
            exe_count += 1
            T2_start = time.time()
            memory = _execute_quantum_program(con_dic['inputID'], con_dic['outputID'], con_dic['num_qubit'], i,
                                                      con_dic['module_name'], shots)
            T2_end = time.time()
            T2 = T2_end - T2_start
            pop_list = []
            for suite in not_finish_list:
                result_dic['input_data'][suite].append(i)
            for o_index in range(len(memory)):
                o = memory[o_index]
                result_dic['output_data'][not_finish_list[o_index]].append(o)
                for uo_index in range(len(com_dic[i])):
                    uo = com_dic[i][uo_index]
                    if o == uo:
                        output_flag[not_finish_list[o_index]][uo_index] = 1
                if _check_same(output_flag[not_finish_list[o_index]], 1):
                    pop_list.append(o_index)
                    # not_finish_list.pop(o_index)
                    shots -= 1
            pop_list = sorted(pop_list, reverse=True)
            for index in pop_list:
                not_finish_list.pop(index)
            if shots == 0:
                break
            if exe_count == BUDGET_I:
                break
    result_dic['T2'] = T2
    return result_dic

def _input_group(valid_input):
    index = [] #unique input index
    index_flag = valid_input[0]
    index.append(0)
    for i in range(1,len(valid_input)):
        if valid_input[i] != index_flag:
            index.append(i)
            index_flag = valid_input[i]
    return index

def _check_full_ps(valid_input, p):
    index = _input_group(valid_input)
    p_sum = 0
    for i in range(len(index)):
        start = index[i]
        if i == len(index) - 1:
            end = len(valid_input)
        else:
            end = index[i+1]
        for j in range(start,end):
            p_sum += p[j]
        if p_sum != 1:
            print("Error: This is not a complete program specification.")
            _end_running()
        else:
            p_sum = 0

def _check_bin(bin_str, n):
    if len(bin_str) != n:
        print("Error: The format of the program specification is wrong.")
        _end_running()
   # print("check bin: "+str(bin_str))
    for i in range(len(bin_str)):
        if bin_str[i] != '0' and bin_str[i] != '1':
            print("Error: The format of the program specification is wrong.")
            _end_running()

def _assign(cate, coverage, con_dic):
    T1_start = time.time()
    resultFolder = os.path.join(con_dic['program_folder'], "result")
    if not os.path.exists(resultFolder):
        os.makedirs(resultFolder)
    if coverage == 'IC':
        input_file_name = '_input_coverage_'
        if cate == 'no':
            result_dic = _input_coverage_no(con_dic)
        elif cate == 'partial':
            result_dic = _input_coverage_partial(con_dic)
        elif cate == 'full':
            result_dic = _input_coverage(con_dic)
    elif coverage == 'OC':
        input_file_name = '_output_coverage_'
        if cate == 'no':
            result_dic = _output_coverage_no(con_dic)
        elif cate == 'partial':
            result_dic = _output_coverage_partial(con_dic)
        elif cate == 'full':
            result_dic = _output_coverage(con_dic)
    elif coverage =='IOC':
        input_file_name = '_input_output_coverage_'
        if cate == 'no':
            result_dic = _input_output_coverage_no(con_dic)
        elif cate == 'partial':
            result_dic = _input_output_coverage_partial(con_dic)
        elif cate == 'full':
            result_dic = _input_output_coverage(con_dic)
    T1_end = time.time()
    T1 = T1_end-T1_start
    input_file = os.path.join(resultFolder, 'INPUTS'+ input_file_name + con_dic['module_name']+'.txt')
    result_file = os.path.join(resultFolder, 'RESULTS'+ input_file_name + con_dic['module_name']+'.txt')
    input_file = open(input_file, "w")
    result_file = open(result_file,'w')
    case_num = 0
    for suite_index in range(len(result_dic['input_data'])):
        for case_index in range(len(result_dic['input_data'][suite_index])):
            try:
                result_file.write('('+result_dic['input_data'][suite_index][case_index]+','+result_dic['output_data'][suite_index][case_index]+')'+' ')
                input_file.write(result_dic['input_data'][suite_index][case_index]+' ')
                case_num += 1
            except:
                continue
        input_file.write('\n')
        result_file.write('\n')
    if cate == 'partial' or cate == 'full':
        ass_file = os.path.join(resultFolder, 'ASSESSMENT'+input_file_name+con_dic['module_name']+'.txt')
        ass_file = open(ass_file,'w')
        if result_dic['WOO']:
            ass_file.write("Unexpected input-output pair: ("+str(list(result_dic['w_input'].keys())[0])+
                           ','+str(list(result_dic['w_input'].values())[0])+')'+'\n')
            ass_file.write("Fail for WOO.")
            result_file.write('\n')
            result_file.write('The total number of test cases is ' + str(case_num) + '.')
            result_file.write('\n')
            result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(result_dic['T2']) + "s.")
        else:
            _judge_ass_result(con_dic['valid_input'], con_dic['valid_output'], result_dic['pvalue'],con_dic['C_LEVEL'], ass_file)
            result_file.write('\n')
            result_file.write('The total number of test cases is ' + str(case_num) + '.')
            result_file.write('\n')
            result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(result_dic['T2']) + "s.")
    print("The total run time is " + "{:.2f}".format(T1) + "s.")
    print('The execution time of the quantum program is ' + "{:.2f}".format(result_dic['T2']) + "s.")

def _quito_run(root_con):
    # global START###
    # START = time.time()###
    #get configuration file

    # Checking the root of configuration file
    if os.path.isfile(root_con) == True:
        config = configparser.ConfigParser(allow_no_value=True)
        try:
            config.read(root_con, encoding='utf-8')
        except:
            print("Error: The configuration has duplicate entries.")
            _end_running()
    else:
        print("Error: Quito cannot find the configuration file.")
        _end_running()

    # Checking the format of configurations and returning the ps category
    ps_category = _check_configuration_file(config)
    if ps_category != 'no' and ps_category != 'partial' and ps_category != 'full':
        print("Error: The format of program specification category is wrong.")
        _end_running()

    # Getting the quantum program
    root = config.get('program','root')
    if not os.path.isfile(root):
        print("Error: Quito cannot find the quantum program file.")
        _end_running()

    program_folder = os.path.dirname(root)  # the root
    program_file = os.path.basename(root)  # the file name
    module_name = program_file.split('.')[0]  # the module name
    sys.path.append(program_folder)

    # root_list = root.split(ROOT)###
    # program_file = root_list[len(root_list)-1]
    # program_folder = root_list[:len(root_list)-1]
    # program_folder = ROOT.join(str(i) for i in program_folder)
    # sys.path.append(program_folder)
    # # print(program_file.split('.')[0])
    # module_name = program_file.split('.')[0]

    #get inputID, outputID and numner of qubits

    # Getting input and output
    inputID_o = config.get('program','inputID').split(',')
    outputID_o = config.get('program','outputID').split(',')
    num_qubit = int(config.get('program','num_qubit'))
    inputID, outputID = _check_inputID_outputID(num_qubit, inputID_o, outputID_o)

    #Getting coverage criterion
    coverage_criterion = config.get('quito_configuration', 'coverage_criterion')
    if coverage_criterion != 'IC' and coverage_criterion != 'OC' and coverage_criterion != 'IOC':
        print("Error: The format of coverage criterion is not right.")
        _end_running()

    # Getting paramaters
    K = 200
    if config.has_option('quito_configuration', 'K') != None:
        try:
            K = int(config.get('quito_configuration', 'K'))
        except:
            print("Error: The data type of 'K' should be an Integer.")
    M = 20
    if config.has_option('quito_configuration', 'M') != None:
        try:
            M = int(config.get('quito_configuration', 'M'))
        except:
            print("Error: The data type of 'M' should be an Integer.")
    C_LEVEL = 0.01
    if config.has_option('quito_configuration', 'confidence_level') != None:
        try:
            C_LEVEL = float(config.get('quito_configuration', 'confidence_level'))
        except:
            print("Error: The data type of 'confidence_level' should be a Float.")
    S_TEST = 'one-sample Wilcoxon signed rank test'

    con_dic = {}
    con_dic['inputID'] = inputID
    con_dic['outputID'] = outputID
    con_dic['num_qubit'] = num_qubit
    con_dic['module_name'] = module_name
    con_dic['program_folder'] = program_folder
    con_dic['K'] = K
    con_dic['M'] = M
    con_dic['C_LEVEL'] = C_LEVEL
    con_dic['S_TEST'] = S_TEST

    if ps_category == 'no':
        if coverage_criterion == 'IC':
            _assign(ps_category, coverage_criterion, con_dic)
            # input_coverage_no(con_dic)
        else:
            con_dic['BUDGET'] = pow(2, len(inputID)) * 10
            if config.has_option('quito_configuration', 'BUDGET') != None:
                try:
                    con_dic['BUDGET'] = int(config.get('quito_configuration', 'BUDGET'))
                except:
                    print("Error: The data type of 'BUDGET' should be an Integer.")
                # check budget
                if con_dic['BUDGET'] < pow(2, len(inputID)):
                    print("Error: Budget is smaller than the number of inputs.")
                    _end_running()
                _assign(ps_category, coverage_criterion, con_dic)
                # if coverage_criterion == 'OC':
                #     output_coverage_no(con_dic)
                # elif coverage_criterion == 'IOC':
                #     input_output_coverage_no(con_dic)
    else:
        # get PS
        valid_input = []
        valid_output = []
        p = []
        ps = config.items('program_specification')
        # print("origin:"+str(ps))

        if _check_unique(ps) == False:
            print("Program specifications not unique")
            _end_running()
        # sort PS according to input and output
        ps.sort(key=lambda x: x[0])
        # print("new:"+str(ps))
        for i in range(len(ps)):
            valid_input_item = ps[i][0][:len(inputID)]
            valid_output_item = ps[i][0][len(inputID) + 1:]
            _check_bin(valid_input_item, len(inputID))
            _check_bin(valid_output_item, len(outputID))
            valid_input.append(valid_input_item)
            valid_output.append(valid_output_item)
            p.append(float(ps[i][1]))

        con_dic['valid_input'] = valid_input
        con_dic['valid_output'] = valid_output
        con_dic['p'] = p

        # check budget
        if coverage_criterion == 'OC' or 'IOC':
            if ps_category == 'full':
                con_dic['BUDGET'] = pow(2, len(set(valid_input))) * 10  # default
                if config.get('quito_configuration', 'BUDGET') != None:
                    con_dic['BUDGET'] = int(config.get('quito_configuration', 'BUDGET'))
                num_valid_input = set(valid_input)
                if con_dic['BUDGET'] < len(num_valid_input):
                    print("Error: Budget is smaller than the number of inputs.")
                    _end_running()
            elif ps_category == 'partial':
                con_dic['BUDGET'] = pow(2, len(inputID)) * 10  # default
                if config.get('quito_configuration', 'BUDGET') != None:
                    con_dic['BUDGET'] = int(config.get('quito_configuration', 'BUDGET'))
                if con_dic['BUDGET'] < pow(2, len(inputID)):
                    print("Error: Budget is smaller than the number of inputs.")
                    _end_running()

        if ps_category == 'full':
            _check_full_ps(valid_input, p)
        # print(outputID)

        _assign(ps_category, coverage_criterion, con_dic)

        # if ps_category == 'full':
        #     if coverage_criterion == 'IC':
        #         _input_coverage(con_dic)
        #     elif coverage_criterion == 'OC':
        #         output_coverage(con_dic)
        #     elif coverage_criterion == 'IOC':
        #         input_output_coverage(con_dic)
        #
        # if ps_category == 'partial':
        #     if coverage_criterion == 'IC':
        #         input_coverage_partial(con_dic)
        #     elif coverage_criterion == 'OC':
        #         output_coverage_partial(con_dic)
        #     elif coverage_criterion == 'IOC':
        #         input_output_coverage_partial(con_dic)
