import pandas as pd
import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm
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
# fig.set_size_inches(5.5, 4.8)
thickness = 7
width = 3
smf = 0.05
intv = 10
per = 10
num = 40

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

Title = 'Cyclic loading(10 loads)_20g_'
### room temp fruit strips
strip1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip1'
strip2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip2'
strip3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip3'
strip4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip4'
# strip_s1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-1_20g'
# strip_s2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-2_20g'
# strip_s3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-3_10g'
# strip_s4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-4_30g'
# strip_s5 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-5_10g'
### strips area
# area_list = [3552762.637, 2927490.149, 3650148.381, 2944460.14]
### standard strips area
# area_list = [3389297.281, 2679108.159, 3082287.14, 3597204.676, 2974816.086]
# area_list = [3082287.14]

### 4 degree fruit strips
# strip1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip1_10g'
# strip2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip2_10g'
# strip3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip3_10g'
# strip4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip4_10g'
# area_list = [4373323.381, 4215275.754, 3970144.4, 3568115.645]
# strip_s1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-1_20g'
# strip_s2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-2_20g'
strip_s3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-3_10g'
# strip_s4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-4_30g'
# strip_s5 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-5_10g'
### standard strips area
# area_list = [3251820.35, 2591324.295, 3610690.651, 3374927.8, 3102991.195]
area_list = [3610690.651]


# target for seuqntial Instron
target = [2,4,6,8,10,10,12,14,16,18,20,20]

Pull = pd.DataFrame(columns= target)
Retract = pd.DataFrame(columns= target)
# group = [strip1, strip2, strip3, strip4]
group = [strip_s3]

# group = [strip_s1, strip_s2, strip_s3, strip_s4, strip_s5]

Strain_4C = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c'
Strain_4C_2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c/Onion3'
Strain_RT_Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_RT_control'

# 20g
Strain_4c_20g = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c_20g'
# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

def pos_p(pull, p_pull):
    for i in range(len(pull) - 5):
        n = np.max(p_pull.strain)
        if (pull.strain[i] > n) & (pull.strain[i + 1] > n) & (pull.strain[i + 2] > n) & \
                (pull.strain[i + 3] > n) & (pull.strain[i + 4] > n):
            index = i
            return index
        elif i == (len(pull) - 6):
            return i

def ten_ind(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            return i

def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i


def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return (len(pull) -1)

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)


# define a function that calculate the stress and strain of the curve
def norm(pull, area_um2):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)
    # strain calculated based on length at the beginning of each pull
    # ori_p = pull.position[0]

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    # pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['stress'] = pull.force_N / area_um2 * 10 ** 9
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # true stress & strain
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p + 4500) )
    # pull['stress'] = pull.force_N/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)

    # when input is strain (new Instron output setting)
    # pull['strain'] = np.log(1 + pull.position)
    # pull['stress'] = pull.force_N/(thickness * width / (pull.position + 1) * 0.001)

def fit_d(pull, num):
    fitting_part = pull[len(pull)-num:].reset_index(drop=True)
    # fit using gram as unit
    # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]
    return z[0]

def fit(pull, percentage, target_load):
    if target_load == 'N':
        fitting_part = pull[int(len(pull) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    else:
        # calculate according to target load
        cutted_data = pull[:tar_ind(pull, target_load)].reset_index(drop=True)
        fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    return z[0]

def fit_p(pull, p_pull, percentage):
    cutted_data = pull[pull.strain < np.max(p_pull.strain)].reset_index(drop=True)
    fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
    # fit using gram as unit
    # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]

    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    return z[0]

# define a function that smooth the curve
def smooth(pull,f):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= f)
    return z

def derv(pull, inter):
    deriv = []
    for i in range(inter):
        deriv.append((pull.stress[i + inter] - pull.stress[i]) / (pull.strain[i + inter] - pull.strain[i]))
    for i in range(inter, len(pull) - inter):
        deriv.append(
            (pull.stress[i + inter] - pull.stress[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
            # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
    for i in range(len(pull) - inter, len(pull)):
        deriv.append((pull.stress[i] - pull.stress[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
    return deriv

def dervuni(y, x, inter):
    deriv = []
    for i in range(inter):
        deriv.append((y[i + inter] - y[i]) / (x[i + inter] - x[i]))
    for i in range(inter, len(y) - inter):
        deriv.append(
            (y[i + inter] - y[i - inter]) / (x[i + inter] - x[i - inter]))
            # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
    for i in range(len(y) - inter, len(y)):
        deriv.append((y[i] - y[i - inter]) / (x[i] - x[i - inter]))
    return deriv

# path = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Nitro_glove_dif_rate_20g/Instron_with_Retract__03_18_2020__15_06_56_SHORT.csv'

#for file in sorted(glob.glob(os.path.join(path, '*SHORT.csv'))):

color = ['C0','C1','C2', 'C3', 'C4', 'C5','C6','C7','C8','C9','C10','C11']
n = 0
cc = 0
for folder in group:
    area = area_list[cc]
    n = 0
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):

        if n%2 == 0:
            print(n)
            whole_df = pd.read_table(file,
                                 delimiter=',',
                                 header=None,
                                 names=range(0, 2),
                                 )
            whole_df.columns = ['position', 'load']


            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)

            ori_p = ten_p(first_pull)
            n += 1
            continue

        if n%2 == 1:
            print(n)
            whole_df = pd.read_table(file,
                                 delimiter=',',
                                 header=None,
                                 names=range(0, 2),
                                 )
            whole_df.columns = ['position', 'load']


            seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
            seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
            eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(drop=True)
            ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(drop=True)
            tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(drop=True)
            eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(drop=True)
            twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)


            # pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull, seventh_pull, eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            # retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract, seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]
            # pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull, eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            # retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

            pre_p = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull, seventh_pull]
            pre_r = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract, seventh_retract]
            rate_p = [eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            rate_r = [eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

            # pull = [first_pull, second_pull]
            # retract = [first_retract, second_retract]

            ### normalize curve
            for i in range(len(pre_p)):
                norm(pre_p[i], area)
                norm(pre_r[i], area)
                # print((ten_p(pre_p[i]) - ori_p)/5000)
            for i in range(len(rate_p)):
                norm(rate_p[i], area)
                norm(rate_r[i], area)
                # print(ten_p(rate_p[i])/5000)

            for i in range(len(pre_r)):
                pre_r[i].drop(pre_r[i].loc[purify_p(pre_r[i], 10): len(pre_r[i])].index, inplace=True)
            for i in range(len(rate_r)):
                rate_r[i].drop(rate_r[i].loc[purify_p(rate_r[i], 10): len(rate_r[i])].index, inplace=True)



            sm1pull = pd.DataFrame(data = smooth(first_pull,smf),columns=['strain','stress'])
            sm2pull = pd.DataFrame(data = smooth(second_pull,smf),columns=['strain','stress'])
            sm3pull = pd.DataFrame(data = smooth(third_pull,smf),columns=['strain','stress'])
            sm4pull = pd.DataFrame(data = smooth(fourth_pull,smf),columns=['strain','stress'])
            sm5pull = pd.DataFrame(data = smooth(fifth_pull,smf),columns=['strain','stress'])
            sm6pull = pd.DataFrame(data = smooth(sixth_pull,smf),columns=['strain','stress'])
            sm7pull = pd.DataFrame(data = smooth(seventh_pull,smf),columns=['strain','stress'])
            sm8pull = pd.DataFrame(data = smooth(eighth_pull,smf),columns=['strain','stress'])
            sm9pull = pd.DataFrame(data = smooth(ninth_pull,smf),columns=['strain','stress'])
            sm10pull = pd.DataFrame(data = smooth(tenth_pull,smf),columns=['strain','stress'])
            sm11pull = pd.DataFrame(data = smooth(eleventh_pull,smf),columns=['strain','stress'])
            sm12pull = pd.DataFrame(data = smooth(twlveth_pull,smf),columns=['strain','stress'])

            sm1retract = pd.DataFrame(data=smooth(first_retract, smf), columns=['strain', 'stress'])
            sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
            sm3retract = pd.DataFrame(data=smooth(third_retract, smf), columns=['strain', 'stress'])
            sm4retract = pd.DataFrame(data=smooth(fourth_retract, smf), columns=['strain', 'stress'])
            sm5retract = pd.DataFrame(data=smooth(fifth_retract, smf), columns=['strain', 'stress'])
            sm6retract = pd.DataFrame(data = smooth(sixth_retract,smf),columns=['strain','stress'])
            sm7retract = pd.DataFrame(data=smooth(seventh_retract, smf), columns=['strain', 'stress'])
            sm8retract = pd.DataFrame(data=smooth(eighth_retract, smf), columns=['strain', 'stress'])
            sm9retract = pd.DataFrame(data=smooth(ninth_retract, smf), columns=['strain', 'stress'])
            sm10retract = pd.DataFrame(data=smooth(tenth_retract, smf), columns=['strain', 'stress'])
            sm11retract = pd.DataFrame(data=smooth(eleventh_retract, smf), columns=['strain', 'stress'])
            sm12retract = pd.DataFrame(data=smooth(twlveth_retract, smf), columns=['strain', 'stress'])

            smpre_p = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull, sm7pull]
            smpre_r = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract,sm7retract]
            smrate_p = [sm8pull, sm9pull, sm10pull, sm11pull, sm12pull]
            smrate_r = [sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]


            # calculate the original position shift of the curve
            # step = [0]
            # for i in range(len(pull)):
            #     if i > 0:
            #         step.append((ten_p(pull[i]) - ten_p(pull[0])) / 5000)
            map = 'Paired'
            a = plt.subplot(111)
            for j in range(len(pre_p)):
                a.plot(pre_p[j].strain *100, pre_p[j].stress, color='grey', alpha = 0.4)
                a.plot(pre_r[j].strain *100, pre_r[j].stress, color='grey', linestyle = '--', alpha = 0.4)
            for j in range(len(rate_p)):
                a.plot(rate_p[j].strain *100, rate_p[j].stress, color= color[j])
                a.plot(rate_r[j].strain *100, rate_r[j].stress, color=color[j], linestyle='--')

            ##### strain-modulus plot
            # b = plt.subplot(312)
            # for j in range(len(pre_p)):
            #     b.plot(smpre_p[j].strain *100, derv(smpre_p[j], intv), color='grey', alpha = 0.4)
            #     b.plot(smpre_r[j].strain *100, derv(smpre_r[j],intv),  color='grey', linestyle = '--', alpha = 0.4)
            #
            # for j in range(len(rate_p)):
            #     b.plot(smrate_p[j].strain *100, derv(smrate_p[j],intv), color= color[j])
            #     b.plot(smrate_r[j].strain *100, derv(smrate_r[j],intv), color=color[j], linestyle='--')

            # b.set(ylim = [-5,400])


            ##### stress-modulus plot
            # c = plt.subplot(212)
            # for j in range(len(pre_p)):
            #     c.plot(smpre_p[j].stress, derv(smpre_p[j], intv), color='grey', alpha = 0.4)
            #     c.plot(smpre_r[j].stress, derv(smpre_r[j],intv),  color='grey', linestyle = '--', alpha = 0.4)
            #
            # for j in range(len(rate_p)):
            #     c.plot(smrate_p[j].stress, derv(smrate_p[j],intv), color= color[j])
            #     c.plot(smrate_r[j].stress, derv(smrate_r[j],intv), color=color[j], linestyle='--')


            # c.set(ylim = [-5,400])
            # loop for seuqntial Instron
            # for i in range(len(pull)):
            #     retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace = True)
            #     retract[i].drop(retract[i].loc[0:ten_ind(retract[i])].index, inplace = True)
            #
            #     retract[i].reset_index(inplace = True)
            #     rmindex = []
            #     purify_revrs(retract[i])
            #     retract[i].drop(retract[i].loc[rmindex].index, inplace = True)
            #
            #     pull[i].drop(pull[i].loc[0:ten_ind(pull[i])].index, inplace = True)
            #     retract[i].reset_index(inplace = True)
            #     pull[i].reset_index(inplace=True)
            # loop for repetitive pulling
            # for i in range(len(pull)):
            #     retract[i].drop(retract[i].loc[purify_p(retract[i], 20): len(retract[i])].index, inplace=True)

                # print(retract[i])
            # first_retract.drop(first_retract.index[[purify_p(first_retract,4),len(first_retract)-1]], inplace = True)

            ### normalize curve
            # for i in range(len(pull)):
            #     norm(pull[i])
            #     norm(retract[i])



            # sm1pull = pd.DataFrame(data = smooth(first_pull,smf),columns=['strain','stress'])
            # sm2pull = pd.DataFrame(data = smooth(second_pull,smf),columns=['strain','stress'])
            # sm3pull = pd.DataFrame(data = smooth(third_pull,smf),columns=['strain','stress'])
            # sm4pull = pd.DataFrame(data = smooth(fourth_pull,smf),columns=['strain','stress'])
            # sm5pull = pd.DataFrame(data = smooth(fifth_pull,smf),columns=['strain','stress'])
            # # sm6pull = pd.DataFrame(data = smooth(sixth_pull,smf),columns=['strain','stress'])
            # sm7pull = pd.DataFrame(data = smooth(seventh_pull,smf),columns=['strain','stress'])
            # sm8pull = pd.DataFrame(data = smooth(eighth_pull,smf),columns=['strain','stress'])
            # sm9pull = pd.DataFrame(data = smooth(ninth_pull,smf),columns=['strain','stress'])
            # sm10pull = pd.DataFrame(data = smooth(tenth_pull,smf),columns=['strain','stress'])
            # sm11pull = pd.DataFrame(data = smooth(eleventh_pull,smf),columns=['strain','stress'])
            # sm12pull = pd.DataFrame(data = smooth(twlveth_pull,smf),columns=['strain','stress'])
            #
            # sm1retract = pd.DataFrame(data=smooth(first_retract,smf), columns=['strain','stress'])
            # sm2retract = pd.DataFrame(data=smooth(second_retract,smf), columns=['strain','stress'])
            # sm3retract = pd.DataFrame(data = smooth(third_retract,smf),columns=['strain','stress'])
            # sm4retract = pd.DataFrame(data = smooth(fourth_retract,smf),columns=['strain','stress'])
            # sm5retract = pd.DataFrame(data = smooth(fifth_retract,smf),columns=['strain','stress'])
            # # sm6retract = pd.DataFrame(data = smooth(sixth_retract,smf),columns=['strain','stress'])
            # sm7retract = pd.DataFrame(data=smooth(seventh_retract,smf), columns=['strain','stress'])
            # sm8retract = pd.DataFrame(data=smooth(eighth_retract,smf), columns=['strain','stress'])
            # sm9retract = pd.DataFrame(data = smooth(ninth_retract,smf),columns=['strain','stress'])
            # sm10retract = pd.DataFrame(data = smooth(tenth_retract,smf),columns=['strain','stress'])
            # sm11retract = pd.DataFrame(data = smooth(eleventh_retract,smf),columns=['strain','stress'])
            # sm12retract = pd.DataFrame(data = smooth(twlveth_retract,smf),columns=['strain','stress'])
            #
            # smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm7pull, sm8pull, sm9pull, sm10pull, sm11pull, sm12pull]
            # smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm7retract, sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]
            # smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull]
            # smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract]

            # ela_target = ['N', 2, 4, 6, 8, 'N', 10, 12, 14, 16, 18, 'N']
            # mo_pull = []
            # mo_re = []
            # for i in range(len(smpull)):
                # retract fitting according to the percentage of the data
                # mo_re.append(fit(smretract[i], per, 'N'))

                # retract fitting according to the time/amount of data point
                # mo_re.append(fit_d(smretract[i], num))
                # fit according to the target load
                # mo_pull.append(fit(smpull[i], 8, ela_target[i]))
                # fit according to the previous position
                # if i > 0:
                #     mo_pull.append(fit_p(smpull[i], smpull[i-1], per))
                # else:
                #     mo_pull.append(fit_p(smpull[i], smpull[i], per))
            # Pull.loc[len(Pull)] = mo_pull
            # Retract.loc[len(Retract)] = mo_re
            # print(mo_pull)
            # print(mo_re)


            # only plot the log curve
            # plt.suptitle(Title + 'log transformed curve' + '_' + str(n), fontweight = 'bold')
            # cx = plt.subplot(121)
            # for j in range(len(smpull)):
            #     cx.plot(smpull[j].strain, np.log(smpull[j].stress), color=color[j])
            # cx.set(ylabel='Log of stress (log(MPa))', xlabel='Strain')
            # cx.set_title('Log stress-strain curve for loading curve', fontsize=10, fontweight='bold')
            # # cx.set_xlim(left = -6)
            # cx.grid(alpha=0.4, linestyle='--')
            #
            # dx = plt.subplot(122)
            # for j in range(len(smpull)):
            #     dx.plot(smretract[j].strain, np.log(smretract[j].stress), color = color[j], linestyle = '--')
            # dx.set(ylabel='Log of stress (log(MPa))', xlabel='Strain')
            # dx.set_title('Log stress-strain curve for retract curve', fontsize=10, fontweight='bold')
            # dx.grid(alpha=0.4, linestyle='--')

            # choose the pulls for visualization
            # smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull]
            # smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract]
            # smpull = [sm1pull, sm2pull, sm3pull]
            # smretract = [sm1retract, sm2retract, sm3retract]

            # plot the comprehensive figure
            # map = 'Paired'
            # plt.suptitle(Title + '_' + str((n+1)/2), fontweight = 'bold')
            # a = plt.subplot(121)
            # for j in range(len(smpull)):
                 # for target loading fitting plot
                # PL = int(tar_ind(smpull[j], ela_target[j]) * (100 - 8) / 100)
                # PR = tar_ind(smpull[j], ela_target[j])
                 # for target position fitting plot
                # if j > 0:
                #     PR = pos_p(smpull[j], smpull[j-1])
                #
                # else:
                #     PR = pos_p(smpull[j], smpull[j])
                # PL = int(PR * (100-per)/100)
                # a.plot(smpull[j].strain , smpull[j].stress, color = color[j])
                # a.plot(smretract[j].strain, smretract[j].stress, color=color[j], linestyle = '--')

            # a.scatter(smpull[j].strain * 100 , smpull[j].stress, c=derv(smpull[j], intv), cmap = map, s = 0.5)
                 # a.scatter(smretract[j].strain *100, smretract[j].stress, c=derv(smretract[j], intv), cmap = map, s = 4)
                 # a.plot(smpull[j].strain[int(len(smpull[j]) * (100 - per) / 100):] * 100 , smpull[j].fit_e[int(len(smpull[j]) * (100 - per) / 100):], color = 'black')
                 # a.plot(smretract[j].strain[int(len(smretract[j]) * (100 - per) / 100):] * 100, smretract[j].fit_e[int(len(smretract[j]) * (100 - per) / 100):], color = 'black')
                # a.plot(smpull[j].strain[PL:PR] * 100 , smpull[j].fit_e[PL:PR], color = 'Black', linewidth = 2)

                 # a.plot(smretract[j].strain[int(len(smretract[j]) * (100 - per) / 100):] * 100, smretract[j].fit_e[int(len(smretract[j]) * (100 - per) / 100):], color = 'Purple', linewidth = 2)

            a.set(xlim = [-0.4, 7], ylabel = 'Stress(kPa)', xlabel = 'Strain(%)')
            # a.set_title('stress-strain curve for loading curve', fontsize = 10, fontweight = 'bold')
            # a.grid(alpha = 0.4, linestyle = '--')
            # b.grid(alpha = 0.4, linestyle = '--')
            # c.grid(alpha = 0.4, linestyle = '--')
            # setting for gram y axis
            ay = a.twinx()
            ay.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
            B = -1 * 9.8
            T = 13 * 9.8
            ay.set(ylim=[B, T])
            SB = B * 0.001 / area * 10 ** 9
            ST = T * 0.001 / area * 10 ** 9
            a.set(ylim=[SB, ST])
            grey_line = mlines.Line2D([], [], color = 'grey', label = 'Pre-stretch (10mm/min)' )
            C0_line = mlines.Line2D([], [], color = 'C0', label = '10 mm/min' )
            C1_line = mlines.Line2D([], [], color = 'C1', label = '5 mm/min' )
            C2_line = mlines.Line2D([], [], color = 'C2', label = '1 mm/min' )
            C3_line = mlines.Line2D([], [], color = 'C3', label = '0.5 mm/min' )
            C4_line = mlines.Line2D([], [], color = 'C4', label = '0.1 mm/min' )
            a.set(xlabel = 'Strain (%)', ylabel = 'Stress (kPa)')
            a.legend(handles = [grey_line, C0_line, C1_line, C2_line, C3_line, C4_line], loc = 2, fontsize=10, frameon = False)
            # b.set(xlabel = 'Strain (%)', ylabel = 'Modulus (MPa)')
            # c.set(xlabel = 'Stress (kPa)', ylabel = 'Modulus (kPa)')


        # b = plt.subplot(122)
            # for j in range(len(smpull)):
                # b.plot(smretract[j].strain , smretract[j].stress, color = color[j], linestyle = '--')
                # b.scatter(smretract[j].strain *100, smretract[j].stress, c=derv(smretract[j], intv), cmap = map, s =2)
                # b.plot(smretract[j].strain[int(len(smretract[j]) * (100 - per) / 100):] * 100, smretract[j].fit_e[int(len(smretract[j]) * (100 - per) / 100):], color = 'black')
                # b.plot(smretract[j].strain[int(len(smretract[j]) * (100 - per) / 100):] * 100, smretract[j].fit_e[int(len(smretract[j]) * (100 - per) / 100):], color = 'black')
                # b.plot(smretract[j].strain[len(smretract)-num:] * 100 , smretract[j].fit_e[len(smretract)-num:], color = 'Black', linewidth = 2)

            # b.set(ylim = [-1,10], ylabel = 'Stress(MPa)', xlabel = 'Strain(%)')
            # b.set_title('stress-strain curve for retract curve', fontsize = 10, fontweight = 'bold')
            # b.grid(alpha = 0.4, linestyle = '--')
            # setting for gram y axis
            # by = b.twinx()
            # by.set_ylabel('Load (grams)')
            # B = -1
            # T = 21
            # by.set(ylim=[B, T])
            # SB = B * 0.0098 / (thickness * width * 0.001)
            # ST = T * 0.0098 / (thickness * width * 0.001)
            # b.set(ylim=[SB, ST])

            # a0 = plt.subplot(312)
            # for j in range(len(smpull)):
                # a0.plot(smpull[j].strain *100, smpull[j].stress, color = color[j])
                # a0.scatter(smpull[j].strain *100, derv(smpull[j], intv), c=derv(smpull[j], intv), cmap = map, s =2)
                # a0.scatter(smpull[j].strain *100, smpull[j].stress, c=derv(smpull[j], intv), cmap = map, s =2)

            # a0.set(ylabel = 'Stress (MPa)', xlabel = 'Strain')
            # a0.set_title('stress-strain curve for loading curve (log)', fontsize = 10, fontweight = 'bold')
            # a0.grid(alpha = 0.4, linestyle = '--')
            # a0.set_yscale('log')

            # b0 = plt.subplot(426)
            #
            # for j in range(len(smpull)):
                # b0.plot(smretract[j].stress, derv(smretract[j], intv), color = color[j], linestyle = '--')
                # b0.scatter(smretract[j].strain *100, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map, s =2)
                # b0.scatter(smretract[j].strain *100, smretract[j].stress, c=derv(smretract[j], intv), cmap = map, s =2)

            # b0.set(ylabel = 'Stress (MPa)', xlabel = 'Strain(%)')
            # b0.set_title('stress-strain curve for retract curve (log)', fontsize = 10, fontweight = 'bold')
            # b0.grid(alpha = 0.4, linestyle = '--')
            # b0.set_yscale('log')

            # ax = plt.subplot(425)
            # for j in range(len(pull)):
            #     # ax.plot(smpull[i].strain *100, smpull[i].stress, color = color[i])
            #     ax.scatter(smpull[j].stress, derv(smpull[j], intv), c=derv(smpull[j], intv), cmap = map, s =2)
            #
            # ax.set(ylabel = 'Derivative (MPa)', xlabel = 'Stress(MPa)')
            # ax.set_title('∆σ/∆ε over σ for loading curve', fontsize = 10, fontweight = 'bold')
            # ax.grid(alpha = 0.4, linestyle = '--')

            # bx = plt.subplot(426)
            # #
            # for j in range(len(smpull)):
            #     # bx.plot(smretract[j].stress, derv(smretract[j], intv), color = color[j], linestyle = '--')
            #     bx.scatter(smretract[j].stress, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map, s =2)
            #
            # bx.set(ylabel = 'Derivative (MPa)', xlabel = 'Stress(MPa)')
            # bx.set_title('∆σ/∆ε over σ for retract curve', fontsize = 10, fontweight = 'bold')
            # bx.grid(alpha = 0.4, linestyle = '--')

            # cx = plt.subplot(423)
            # for j in range(len(smpull)):
            #     cx.plot(smpull[j].strain , smpull[j].stress, color = color[j])
            #     cx.scatter(smpull[j].strain * 100, derv(smpull[j], intv), c=derv(smpull[j], intv), cmap = map,s = 0.5)
                # cx.scatter(smretract[j].strain , derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map,s =2)
                # cx.scatter(smpull[j].strain , derv(smpull[j],intv), c=derv(smpull[j], intv), s =2)
            #
            # cx.set(ylim = [-5,130], ylabel = 'Derivative (MPa)', xlabel = 'Strain(%)')
            # cx.set_title('derivative (loading curve)', fontsize = 10, fontweight = 'bold')
            #
            # cx.grid(alpha = 0.4, linestyle = '--')
            # cx.set_yscale('log')


            # dx = plt.subplot(424)
            # for j in range(len(smpull)):
                # dx.plot(smretract[j].strain , np.log(smretract[j].stress), color = color[j], linestyle = '--')
                # dx.scatter(smretract[j].strain *100, smretract[j].stress, c=derv(smretract[j], intv), cmap = map, s =2)
                # dx.scatter(smretract[j].strain * 100, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map,s = 0.5)
                # dx.scatter(smretract[j].strain , np.log(derv(smretract[j], intv)), c=derv(smretract[j], intv), s =2)

            # dx.set(ylim = [-5,500], ylabel = 'Derivative(MPa)', xlabel = 'Strain(%)')
            # dx.set_title('derivative (retract curve)', fontsize = 10, fontweight = 'bold')
            # dx.grid(alpha = 0.4, linestyle = '--')
            # dx.set_yscale('log')

            # ex = plt.subplot(313)
            # for j in range(len(smpull)):
            #     ex.plot(smretract[j].strain , np.log(smretract[j].stress), color = color[j], linestyle = '--')
            #     ex.scatter(smpull[j].strain * 100, dervuni(np.log(smpull[j].stress), smpull[j].strain, intv), c=derv(smpull[j], intv), cmap=map, s=2)
                # ex.scatter(smretract[j].strain * 100, dervuni(np.log(smretract[j].stress), smretract[j].strain, intv), c=derv(smretract[j], intv), cmap=map, s=2)
                # ex.scatter(smretract[j].strain * 100, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map,s =2)

                # ex.scatter(smretract[j].strain , np.log(derv(smretract[j], intv)), c=derv(smretract[j], intv), s =2)

            # ex.set(ylim = [-5, 300], ylabel='Derivative (MPa)', xlabel='Strain(%)')
            # ex.set_title('Log derivative (loading curve)', fontsize=10, fontweight='bold')
            # ex.grid(alpha=0.4, linestyle='--')
            # ex.set_yscale('log')


            # fx = plt.subplot(428)
            # for j in range(len(smpull)):
                # fx.plot(smretract[j].strain , np.log(smretract[j].stress), color = color[j], linestyle = '--')
                # fx.scatter(smpull[j].strain * 100, dervuni(np.log(smpull[j].stress), smpull[j].strain, intv), c=derv(smpull[j], intv), cmap=map, s=2)
                # fx.scatter(smretract[j].strain * 100, dervuni(np.log(smretract[j].stress), smretract[j].strain, intv), c=derv(smretract[j], intv), cmap=map, s=2)
                # fx.scatter(smretract[j].strain * 100, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map,s =2)

                # fx.scatter(smretract[j].strain , np.log(derv(smretract[j], intv)), c=derv(smretract[j], intv), s =2)

            # fx.set(ylim = [-5, 300], ylabel='Derivative (MPa)', xlabel='Strain(%)')
            # fx.set_title('Log derivative (retract curve)', fontsize=10, fontweight='bold')
            # fx.grid(alpha=0.4, linestyle='--')
            n += 1
            cc += 1
        # if n == 4:
        #     break
        #     plt.show()
        #     break

            plt.gcf().subplots_adjust(bottom=0.15, right= 0.815)
            plt.savefig(''f'{output}/Strain_rate_dependency(pectin_4C).pdf', transparent=True)
            plt.show()

# mean_pull = []
# mean_retract = []
# for i in range(len(target)):
#     mean_pull.append(np.mean(Pull.iloc[:, i]))
#     mean_retract.append(np.mean(Retract.iloc[:, i]))
# print(Pull, Retract)
#
# jx = plt.subplot(211)
# jx.bar(target, mean_pull)
# px = plt.subplot(212)
# px.bar(target, mean_retract)

# plt.show()
# plt.bar(target, Retract)

# Pull.to_csv('/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/Pull_P10.csv')
# Retract.to_csv('/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/Retract_P10.csv')
# ax = plt.subplot(111)
# blue_line = mlines.Line2D([], [], color = 'blue', label = '1˚cycle' )
# green_line = mlines.Line2D([], [], color = 'green', label = '2˚cycle' )
# orange_line = mlines.Line2D([], [], color = 'orange', label = '3˚cycle' )
# black_line = mlines.Line2D([], [], color = 'black', label = '4˚cycle' )
# ax.set(xlabel = 'strain', ylabel = 'Force (N)', title = 'Nitro-glove cyclic pulling_10mm/min')
# ax.legend(handles = [blue_line, green_line, orange_line, black_line])




