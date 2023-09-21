# script: Study progression of elasticity and plasticity with compliance method
# output value of T & E & P compliance and their incremental compliance;
# plot: line-plot for exact value, stack bar chart for incremental value

# import necessary package
import pandas as pd
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from math import e
from numpy import trapz
import statsmodels.api as sm
import copy
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.use('macosx')
fig = mpl.pyplot.gcf()
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

# User input
smf = 0.04
folder = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/10mmpermin_9.22.20'
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.17.20'
folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion1'
Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
Onion4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion4'
Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'
# Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/Con_pH7.0'
# PEG = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/40percentPEG8k_pH7.0'
example = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load'

# de-waxed onion peel
Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Chloroform_methanol_treated_onion/Con_freshedly_peeled_washed'
Chlo = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Chloroform_methanol_treated_onion/Chloroform_ONtreated'
Chlo_met = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Chloroform_methanol_treated_onion/Chloroform_mechanol_1vs1_ONtreated'

# group = [folder, folder2, folder3, folder4]
# group = [folder2, Onion3, Onion5]
group = [example, folder2, Onion3, Onion5, Onion6]
# group = [Con, Chlo, Chlo_met]
# group = [folder2]

# total = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Total_3onion'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/seq_Onion_10'
# folder5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Onion_6.23.20_9'
# group = [folder1, folder2, folder3, folder4, folder5]

# 8,12,16,20,24 sequence
seq1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/Onion1'
seq2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/onion2_10'

name = 'Onion4'
# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output = ''f'{folder1}/summary1.csv'

thickness = 7
width = 3
# linear fitting percantage (%)
fitpc = 5

# for load sequence of 4-8-12-16-20
loadx = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
# for load sequence of 8-12-16-20-24
# loadx = [8, 12, 16, 20, 24]


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i


# define a function that find the index at certain loads
def ten_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        for i in range(len(pull) - 5):
            if (pull.load[i] > load) & (pull.load[i + 1] > load) & (pull.load[i + 2] > load) & \
                    (pull.load[i + 3] > load) & (pull.load[i + 4] > load):
                index = i
                return index


def pos_p(pull, p_pull):
    for i in range(len(pull) - 5):
        n = np.max(p_pull.strain)
        if (pull.strain[i] > n) & (pull.strain[i + 1] > n) & (pull.strain[i + 2] > n) & \
                (pull.strain[i + 3] > n) & (pull.strain[i + 4] > n):
            index = i
            return index
        elif i == (len(pull) - 6):
            return i

# define a function that find the position where peels are under tension ( constant > 0.1g )
def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

# define a function that smooth the curve
def smooth(pull,f):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= f)
    return z

# define a function that calculate the stress and strain of the curve
def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)
    # strain calculated based on length at the beginning of each pull
    # ori_p = ten_p(pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # true stress & strain
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p + 4500) )
    # pull['stress'] = pull.force_N/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)

def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return len(pull) -1

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)

# extract the fitting part of the data and calculate the compliance according to the fitting part
# (option: fit using load/stress)
def fit(pull, percentage, target_load):
    if target_load == 'N':
        fitting_part = pull[int(len(pull) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    else:
        cutted_data = pull[:ten_ind(pull, target_load)].reset_index(drop=True)
        fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    modulus = z[0]
    return modulus

def fit_p(pull, p_pull, percentage):
    cutted_data = pull[pull.strain < np.max(p_pull.strain)].reset_index(drop=True)
    fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
    # fit using gram as unit
    # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]

    # Xm = np.median(fitting_part.position)
    # modulus = 100 / (4500 + Xm) / z[0]
    modulus = z[0]
    return modulus

def fit_expo(pull, weight):
    ## fit accroding to original data
    z = np.polyfit(pull.strain, np.log(pull.stress), 1, w = weight)
    pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain)
    ## fit accroding to subtraction of extended linear value
    # z = np.polyfit(pull.strain, np.log(pull.stress - (pull.strain * load_lnr_0.loc[f,loadx[markj-1]] + load_lnr_1.loc[f,loadx[markj-1]])), 1, w = weight)
    # pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain) + pull.strain * load_lnr_0.loc[f,loadx[markj-1]] + load_lnr_1.loc[f,loadx[markj-1]]

    equation = [z[0], z[1]]
    return equation

def fit_expo_l(pull, weight):
    ## fit accroding to original data
    z = np.polyfit(pull.strain, np.log(pull.stress), 1, w = weight)
    pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain)
    ## fit accroding to subtraction of extended linear value
    # z = np.polyfit(pull.strain, np.log(pull.stress - (pull.strain * load_lnr_0.loc[f,loadx[markj-1]] + load_lnr_1.loc[f,loadx[markj-1]])), 1, w = weight)
    # pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain) + pull.strain * load_lnr_0.loc[f,loadx[markj-1]] + load_lnr_1.loc[f,loadx[markj-1]]

    equation = [z[0], z[1]]
    return equation


def fit_linear(pull, weight):
    z = np.polyfit(pull.strain, pull.stress, 1, w = weight)
    pull['fit_e'] = z[0]* pull.strain + z[1]

    equation = [z[0], z[1]]
    return equation

def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i

def Rsquare(curve):
    mean = curve.stress.mean()
    sum_total = 0
    sum_residue = 0
    r = 0
    for i in range(len(curve)):
        sum_total += (curve.stress[i] - mean)**2
        sum_residue += (curve.stress[i] - curve.fit_e[i])**2
    r = 1 - sum_residue/sum_total
    return r



# extract data from each pull
# Use if multiple files need analysis
g = 0
mark = '8'
total = pd.DataFrame(columns=loadx)
ela = pd.DataFrame(columns=loadx)
retra = pd.DataFrame(columns=loadx)
# for folder in group:
#
#     # create data frame for results from each property
#     t_modulus = pd.DataFrame(columns=loadx)
#     e_modulus = pd.DataFrame(columns=loadx)
#     r_modulus = pd.DataFrame(columns=loadx)
#     # t_inc_comp = pd.DataFrame(columns=loadx)
#     # e_inc_comp = pd.DataFrame(columns=loadx)
#     # r_inc_comp = pd.DataFrame(columns=loadx)
#
#     if g == 0:
#         mark = '8'
#     elif g == 1:
#         mark = 's'
#     elif g == 2:
#         mark = 'p'
#     elif g == 3:
#         mark = 'P'
#
n = 0
color = ['C0','C1','C2', 'C3', 'C4', 'C5','C6','C8','C9','C10','C11', 'C12']
retract_exp_0 = pd.DataFrame(columns=loadx)
retract_exp_1 = pd.DataFrame(columns=loadx)
load_exp_0 = pd.DataFrame(columns=loadx)
load_exp_1 = pd.DataFrame(columns=loadx)
load_lnr_0 = pd.DataFrame(columns=loadx)
load_lnr_1 = pd.DataFrame(columns=loadx)
f = 0
markj = 0
for folder in group:
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        if n % 2 == 0:
            # print(n)
            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
                float).reset_index(drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
                float).reset_index(drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
                float).reset_index(drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
                float).reset_index(drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
                drop=True)

            ori_p = ten_p(first_pull)
            n += 1
            continue

        if n % 2 == 1:
            # print(n)
            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
                float).reset_index(drop=True)
            seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
                float).reset_index(drop=True)
            ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
                float).reset_index(drop=True)
            eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
                float).reset_index(drop=True)
            twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
                drop=True)

            n += 1

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull, eighth_pull,
                    ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,
                       seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

        # calculate the stress and strain of the curve
        for i in range(len(pull)):
            norm(pull[i])
            norm(retract[i])

        target = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20]
        ela_target = ['N', 2, 4, 6, 8, 10, 12, 14, 16, 18, 'N']
        ela1_epon = first_pull
        ela2_epon = second_pull
        ela3_epon = third_pull
        ela4_epon = fourth_pull
        ela5_epon = fifth_pull
        ela7_epon = seventh_pull
        ela8_epon = eighth_pull
        ela9_epon = ninth_pull
        ela10_epon = tenth_pull
        ela11_epon = eleventh_pull
        ela12_epon = twlveth_pull

        ela1_lnr = first_pull
        ela2_lnr = second_pull
        ela3_lnr = third_pull
        ela4_lnr = fourth_pull
        ela5_lnr = fifth_pull
        ela7_lnr = seventh_pull
        ela8_lnr = eighth_pull
        ela9_lnr = ninth_pull
        ela10_lnr = tenth_pull
        ela11_lnr = eleventh_pull
        ela12_lnr = twlveth_pull

        re_fit1 = first_retract
        re_fit2 = second_retract
        re_fit3 = third_retract
        re_fit4 = fourth_retract
        re_fit5 = fifth_retract
        re_fit7 = seventh_retract
        re_fit8 = eighth_retract
        re_fit9 = ninth_retract
        re_fit10 = tenth_retract
        re_fit11 = eleventh_retract
        re_fit12 = twlveth_retract

        orip1 = copy.deepcopy(first_pull)
        orip2 = copy.deepcopy(second_pull)
        orip3 = copy.deepcopy(third_pull)
        orip4 = copy.deepcopy(fourth_pull)
        orip5 = copy.deepcopy(fifth_pull)
        orip7 = copy.deepcopy(seventh_pull)
        orip8 = copy.deepcopy(eighth_pull)
        orip9 = copy.deepcopy(ninth_pull)
        orip10 = copy.deepcopy(tenth_pull)
        orip11 = copy.deepcopy(eleventh_pull)
        orip12 = copy.deepcopy(twlveth_pull)

        orir1 = copy.deepcopy(first_retract)
        orir2 = copy.deepcopy(second_retract)
        orir3 = copy.deepcopy(third_retract)
        orir4 = copy.deepcopy(fourth_retract)
        orir5 = copy.deepcopy(fifth_retract)
        orir7 = copy.deepcopy(seventh_retract)
        orir8 = copy.deepcopy(eighth_retract)
        orir9 = copy.deepcopy(ninth_retract)
        orir10 = copy.deepcopy(tenth_retract)
        orir11 = copy.deepcopy(eleventh_retract)
        orir12 = copy.deepcopy(twlveth_retract)


        elaexpo = [ela1_epon, ela2_epon,ela3_epon, ela4_epon, ela5_epon, ela7_epon, ela8_epon, ela9_epon, ela10_epon, ela11_epon, ela12_epon]
        elalnr = [ela1_lnr, ela2_lnr, ela3_lnr, ela4_lnr, ela5_lnr, ela7_lnr, ela8_lnr, ela9_lnr, ela10_lnr, ela11_lnr, ela12_lnr]
        retract_fit = [re_fit1, re_fit2, re_fit3, re_fit4, re_fit5, re_fit7, re_fit8, re_fit9, re_fit10, re_fit11, re_fit12]
        orip = [orip1, orip2, orip3, orip4, orip5, orip7, orip8, orip9, orip10, orip11, orip12]
        orir = [orir1, orir2, orir3, orir4, orir5, orir7, orir8, orir9, orir10, orir11, orir12]

    # loop for seuqntial Instron
        for i in range(len(pull)):
            retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)
            retract[i].reset_index(inplace=True)
            retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
            pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace = True)
            pull[i].reset_index(inplace=True)


            #extract elastic part of the curve
            pull[i].drop(pull[i].loc[ten_ind(pull[i],ela_target[i]):].index, inplace = True)
            pull[i].reset_index(inplace=True)
            if i > 1:
                elalnr[i] = pull[i].loc[int(len(pull[i])*0.1):int(len(pull[i])*0.35)]
            elaexpo[i] = pull[i].loc[int(len(pull[i])*0.35):int(len(pull[i])*0.88)]


            retract[i].reset_index(inplace=True)
            rmindex = []
            # purify_revrs(retract[i])
            # retract[i].drop(retract[i].loc[rmindex].index, inplace=True)
            retract_fit[i] = retract[i].loc[int(len(retract[i])*0.2):int(len(retract[i])*0.95)]

            retract_fit[i] = retract_fit[i].reset_index(drop=True)
            elaexpo[i] = elaexpo[i].reset_index(drop=True)
            elalnr[i] = elalnr[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)

        sm1retract = pd.DataFrame(data=smooth(first_retract,smf), columns=['strain','stress'])
        sm2retract = pd.DataFrame(data=smooth(second_retract,smf), columns=['strain','stress'])
        sm3retract = pd.DataFrame(data = smooth(third_retract,smf),columns=['strain','stress'])
        sm4retract = pd.DataFrame(data = smooth(fourth_retract,smf),columns=['strain','stress'])
        sm5retract = pd.DataFrame(data = smooth(fifth_retract,smf),columns=['strain','stress'])
        sm7retract = pd.DataFrame(data = smooth(seventh_retract,smf),columns=['strain','stress'])
        sm8retract = pd.DataFrame(data = smooth(eighth_retract,smf),columns=['strain','stress'])
        sm9retract = pd.DataFrame(data = smooth(ninth_retract,smf),columns=['strain','stress'])
        sm10retract = pd.DataFrame(data = smooth(tenth_retract,smf),columns=['strain','stress'])
        sm11retract = pd.DataFrame(data = smooth(eleventh_retract,smf),columns=['strain','stress'])
        sm12retract = pd.DataFrame(data = smooth(twlveth_retract,smf),columns=['strain','stress'])


        smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm7retract, sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]

    ####### plot the stress-train curve and fitting curve

        ax = plt.subplot(221)
        bx = plt.subplot(223)
        cx = plt.subplot(222)
        dx = plt.subplot(426)
        dx2 = plt.subplot(428)
        for j in range(len(retract_fit)-1):
            weight_re = []
            for i in range(len(retract_fit[j])):
                if i <= int(len(retract_fit[j]) * 0.3):
                    weight_re.append(1)
                elif int(len(retract_fit[j]) * 0.3) < i <= int(len(retract_fit[j]) * 0.95):
                    weight_re.append(1)
                else:
                    weight_re.append(1)

            fitNmnow = fit_expo(retract_fit[j], weight_re)
            ax.plot(retract[j].strain * 100, retract[j].stress, color = color[j])
            # ax.set_yscale('log')
            ax.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)', title = 'Fitting of elastic retract curves')
            retract_fit_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
            ax.legend(handles = [retract_fit_line], loc = 2, fontsize = 8)
            retract_exp_0.loc[f,loadx[j]] = fitNmnow[0]
            retract_exp_1.loc[f,loadx[j]] = fitNmnow[1]
            ax.plot(retract_fit[j].strain * 100, retract_fit[j].fit_e, color = 'grey')
            bx.scatter(retract_fit[j].strain * 100, (retract_fit[j].stress - retract_fit[j].fit_e)/retract_fit[j].stress *100, color = color[j], s=2)
            bx.set(title = 'Residue of exponential fitting', xlabel = 'Strain(%)', ylabel = 'Residue (%)')



            ax.grid(alpha=0.4, linestyle='--')
            bx.grid(alpha=0.4, linestyle='--')
            cx.grid(alpha=0.4, linestyle='--')
            dx.grid(alpha=0.4, linestyle='--')
            dx2.grid(alpha=0.4, linestyle='--')

        #
        # for j in range(len(elaexpo)):
        for j in range(1, len(loadx)+1):
            cx.plot(pull[j].strain * 100, pull[j].stress, color = color[j])
            cx.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)', title = 'Fitting of elastic loading curves')

            weight_re = []
            for i in range(len(elalnr[j])):
                if i < int(len(elalnr[j]) * 0.05):
                    weight_re.append(1)
                else:
                    weight_re.append(1)

            fitl = (fit_linear(elalnr[j], weight_re))
            # print(fitNmnow)
            load_lnr_0.loc[f,loadx[j-1]] = fitl[0]
            load_lnr_1.loc[f,loadx[j-1]] = fitl[1]
            cx.plot(elalnr[j].strain * 100, elalnr[j].fit_e, color = 'black')
            dx.scatter(elalnr[j].strain * 100, (elalnr[j].stress - elalnr[j].fit_e)/elalnr[j].stress *100 ,color = color[j], s=2)
            dx.set(title = 'Residue of linear fitting', xlabel = 'Strain(%)', ylabel = 'Residue (%)')
            load_expo_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
            load_lnr_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting curve')
            cx.legend(handles = [load_lnr_line, load_expo_line], loc = 2, fontsize = 8)

        for j in range(2, len(loadx)+1):
            weight_re = []
            markj = j
            for i in range(len(elaexpo[j])):
                if i > int(len(elaexpo[j]) * 0.9):
                    weight_re.append(1)
                else:
                    weight_re.append(1)

            fitNmnow = (fit_expo_l(elaexpo[j], weight_re))
            load_exp_0.loc[f,loadx[j-1]] = fitNmnow[0]
            load_exp_1.loc[f,loadx[j-1]] = fitNmnow[1]
            print(load_exp_0)
            cx.plot(elaexpo[j].strain * 100, elaexpo[j].fit_e, color = 'grey')
            dx2.scatter(elaexpo[j].strain * 100, (elaexpo[j].stress - elaexpo[j].fit_e)/elaexpo[j].stress *100, color = color[j], s = 2)
            dx2.set(title = 'Residue of exponential fitting', xlabel = 'Strain(%)', ylabel = 'Residue (%)')

        # plt.show()
        fig.set_size_inches(8, 6)

        master = plt.subplot(111)
        lw = 2
        for i in range(len(retract)):
            master.plot(orir[i].strain *100, orir[i].stress, color = color[i],linestyle = '--', linewidth = lw)
            master.plot(orip[i].strain *100, orip[i].stress, color = color[i], linewidth = lw)
        for i in range(len(retract_fit)-1):
            master.plot(retract_fit[i].strain *100, retract_fit[i].fit_e, color = 'grey', linewidth = lw)
        for i in range(1,len(elalnr)):
            master.plot(elalnr[i].strain *100, elalnr[i].fit_e, color = 'black', linewidth = lw)
        for i in range(2,len(elaexpo)):
            master.plot(elaexpo[i].strain *100, elaexpo[i].fit_e, color = 'grey', linewidth = lw)
        master.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)', title = 'Fitting of elastic curves')
        # master.grid(alpha=0.4, linestyle='--')
        # experiment =
        load_expo_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
        load_lnr_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting curve')
        master.legend(handles = [load_lnr_line, load_expo_line], loc = 2, fontsize=10, frameon = False)

        print(f)
        f += 1
        # plt.gcf().subplots_adjust(bottom=0.12, right = 0.8)

        # plt.savefig(''f'{output}/Seq_cyclic_fitting_master.pdf', transparent = True)

        # plt.show()
# print(load_lnr_0)

sa = plt.subplot(211)
sa2 = plt.subplot(212)
# sb = plt.subplot(212)
R_0 = []
R_1 = []
L_l0 = []
L_l1 = []
L_e0 = []
L_e1 = []
R_0ero = []
R_1ero = []
L_l0ero = []
L_l1ero = []
L_e0ero = []
L_e1ero = []
for i in range(len(loadx)):
    R_0.append(np.mean(retract_exp_0.iloc[:, i]))
    # R_1.append(np.mean(retract_exp_1.iloc[:, i]))
    L_l0.append(np.mean(load_lnr_0.iloc[:, i]))
    # L_l1.append(np.mean(load_lnr_1.iloc[:, i]))

    R_0ero.append(np.std(retract_exp_0.iloc[:, i]))
    # R_1ero.append(np.std(retract_exp_1.iloc[:, i]))
    L_l0ero.append(np.std(load_lnr_0.iloc[:, i]))
    # L_l1ero.append(np.std(load_lnr_1.iloc[:, i]))
    if i > 0:
        L_e0.append(np.mean(load_exp_0.iloc[:, i]))
        # L_e1.append(np.std(load_exp_1.iloc[:, i]))
        L_e0ero.append(np.std(load_exp_0.iloc[:, i]))
        # L_e1ero.append(np.std(load_exp_1.iloc[:, i]))
# print(len(load_exp_1))
loadx_lexpo = [4, 6, 8, 10, 12, 14, 16, 18, 20]
xlabel = [0,1,2,3,4,5,6,7,8,9,10]

say2 = sa.twinx()
say2.set_ylabel('Load exponent', rotation=270, va='bottom')
# say2.set_yticks(range(0, 80, 20))
# sa.set_yticks(range(0, 80, 20))
say2.set(ylim = [0, 70])
sa.set(ylim = [0, 50])
# sa2.set_yticks(range(0, 150, 25))
sa2.set(ylim = [0, 150])

# plot as stress
x = [a * 0.098 /3/7/0.01 for a in loadx]
x_expo = [a * 0.098 /3/7/0.01 for a in loadx_lexpo]

ms = 4
lw = 2
mark = '8'
sa2.errorbar(x, R_0, yerr = R_0ero, color= 'gray', marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
# sb.errorbar(loadx, R_1, yerr = R_1ero, color='C7', alpha = 0.6)
sa.errorbar(x, L_l0, yerr = L_l0ero, color='black', marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
# sb.errorbar(loadx, L_l1, yerr = L_l1ero, color='C9', alpha = 0.6)
say2.errorbar(x_expo, L_e0, yerr = L_e0ero, color='gray', marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
# sb.errorbar(loadx_lexpo, L_e1, yerr = L_e1ero, color='C3', alpha = 0.6)
# sa.grid(alpha=0.4, linestyle='--')
# sb.grid(alpha=0.4, linestyle='--')

# ax.set(xlabel='loading (g)', ylabel='Modulus(MPa)', title='Sequential modulus')
t_line = mlines.Line2D([], [], color = 'C7', label = 'Retract expo base')
# green_line = mlines.Line2D([], [], color = 'C7', label = 'Retract expo intercept', linestyle = '--')
orange_line = mlines.Line2D([], [], color = 'C9', label = 'Load linear slope')
# fr_line = mlines.Line2D([], [], color = 'C9', label = 'Load linear intercept', linestyle = '--')
blue_line = mlines.Line2D([], [], color = 'C3', label = 'Load expo base')
# C_line = mlines.Line2D([], [], color = 'C3', label = 'Load expo intercept', linestyle = '--')

sa.legend(handles = [orange_line, blue_line], loc = 1, fontsize=10, frameon = False)
sa.set(xlabel = 'Stress (MPa)')
sa.set_xticks(xlabel)
sa.set_xticklabels(xlabel)

sa2.legend(handles = [t_line], loc = 1, fontsize=10, frameon = False)
sa2.set(xlabel = 'Stress (MPa)')
sa2.set_xticks(xlabel)
sa2.set_xticklabels(xlabel)
# sa.set_yscale('log')
# sb.legend(handles = [green_line, fr_line, C_line], loc = 1, fontsize = 8)
# sb.set(title = 'Secondary parameters', ylim = [-7,5], xlabel = 'Loading (g)')
# sb.set_xticks(loadx)
# sb.set_xticklabels(loadx)
# ax.grid(alpha=0.4, linestyle='--')
plt.gcf().subplots_adjust(bottom=0.12, right = 0.84, top = 0.95, hspace = 0.35)
plt.savefig(''f'{output}/Seq_cyclic_fitting_parameter.pdf', transparent = True)

plt.show()

