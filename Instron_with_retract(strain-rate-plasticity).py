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
fig.set_size_inches(6.4, 4.8)

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

thickness = 7
width = 3
target = 15
inter = 5
Title = ''

# onion1
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/10mm_min'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/1mm_min'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/0.1mm_min'

#onion2
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/10mm_min'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/1mm_min'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/0.1mm_min'

#all onion
#26C
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/0.1'
#4C
R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/10'
R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/1'
R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/0.1'
## RT control
# R01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/2d_inFridge_0.1'

# fresh = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th'
# nextday = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th/the next day test'


### file output directory
# fR10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-10.csv'
# fR1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-1.csv'
# fR01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-0.1.csv'
#4C
# fR10_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-10.csv'
# fR1_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-1.csv'
# fR01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-0.1.csv'

# fnextday = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th/5th.csv'

#1 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/0.1'
#2 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/0.1'
#3 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/0.1'

# 4C experiment
# R10_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/10mm_min'
# R1_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/1mm_min'
# R01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/0.1mm_min'


group = [R10, R1, R01]
# group = [R10, R1, R01, R10_4c, R1_4c, R01_4c]
# group_name = ['R10', 'R1', 'R01', 'R10_4c', 'R1_4c', 'R01_4c']
# output = [fR10, fR1, fR01, fR10_4c, fR1_4c, fR01_4c]
# group = [fresh ,nextday]
# group_name = ['Nextday']
# output = [fnextday]
# group = [R10b, R1b, R01b]
# group = [R10_4c, R1_4c, R01_4c]
# group = [N_Con, N_4C]
# group = [Con, Xylan]
# group = [unabraded, abraded]
# group = [Ed_Con, Ed_A12P1]
# group = [Ed2_con, Ed2_A2P1, Ed2_E6P2, Ed2_G9P2]

# group = [R10c, R1c, R01c]
Slist = np.arange(0, 4.5, 0.05)
R10c = pd.DataFrame(columns=Slist)
R1c = pd.DataFrame(columns=Slist)
R01c = pd.DataFrame(columns=Slist)
R10c2 = pd.DataFrame(columns=Slist)
R1c2 = pd.DataFrame(columns=Slist)
R01c2 = pd.DataFrame(columns=Slist)
R10c_4c = pd.DataFrame(columns=Slist)
R1c_4c = pd.DataFrame(columns=Slist)
R01c_4c = pd.DataFrame(columns=Slist)
R10c_4c2 = pd.DataFrame(columns=Slist)
R1c_4c2 = pd.DataFrame(columns=Slist)
R01c_4c2 = pd.DataFrame(columns=Slist)

# ssFresh = pd.DataFrame(columns=Slist)
# ssNext = pd.DataFrame(columns=Slist)

S_group = [R10c, R1c, R01c]
S_group2 = [R10c2, R1c2, R01c2]

# S_group = [R10c, R1c, R01c, R10c_4c, R1c_4c, R01c_4c]
# S_group = [ssFresh, ssNext]

summary = pd.DataFrame()


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

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

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

# define a function that calculate the stress and strain of the curve
def norm(pull):
    ori_p = ten_p(first_pull)
    ## stress output
    pull['force'] = pull.load * 0.0098
    pull['stress'] = pull.force / (thickness * width * 0.001)
    ## load output
    # pull['stress'] = pull.load
    # for data output as strain (0 as 5mm)
    # pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # for data output as strain (500 as 5mm)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # for data output as strain
    # pull['strain'] = pull.position * 4500/5000

    # pull['T_strain'] = np.log(1 + pull['strain'])
    # pull['T_stress'] = pull.force/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)

# define a function that smooth the curve
def smooth(pull):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain, frac=0.05)
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

def strain(pull, stress_OI):
    for i in range(len(pull)):
        if pull.stress[i] > stress_OI:
            strainT = pull.strain[i]
            return strainT

# path = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Nitro_glove_dif_rate_20g/Instron_with_Retract__03_18_2020__15_06_56_SHORT.csv'

#for file in sorted(glob.glob(os.path.join(path, '*SHORT.csv'))):

color = ['blue', 'orange', 'red']
# color = ['orange', 'red', 'purple', 'blue', 'green', 'black']
# style = ['solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed']
# marker = ['', '', '', 'v', 'v', 'v']
n = 1
g = 0
a = 0.8
# pla_outputlist = pd.DataFrame(columns= group_name)

for folder in group:
    p_defo = pd.DataFrame(columns= ['defo'])

    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
    # if n <= 3:
    #     n += 1
    #     continue
        whole_df = pd.read_table(file,
                             delimiter=',',
                             header=None,
                             names=range(0, 2),
                             )
        whole_df.columns = ['position', 'load']


        first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
        first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).reset_index(drop=True).iloc[::-1]
        second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
        second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:].astype(float).reset_index(drop=True).iloc[::-1]

        # first_retract.drop(first_retract.index[[purify_p(first_retract, target), len(first_retract)-1]], inplace = True)
        # second_retract.drop(second_retract.index[[purify_p(second_retract, target), len(second_retract)-1]], inplace = True)
        ori_p = ten_p(first_pull)
        norm(first_pull)
        # norm(first_retract)
        norm(second_pull)
        # norm(second_retract)


        retract = [first_retract, second_retract]
        pull = [first_pull, second_pull]

        for i in range(len(pull)):
            # retract[i].drop(retract[i].loc[purify_p(retract[i], 10): len(retract[i])].index, inplace=True)
            # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
            # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

            # retract[i].reset_index(inplace=True)
            rmindex = []
            # purify_revrs(retract[i])
            # retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

            # retract[i] = retract[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)

        sm1pull = pd.DataFrame(data=smooth(first_pull), columns=['strain', 'stress'])
        # sm1retract = pd.DataFrame(data=smooth(first_retract), columns=['strain', 'stress'])
        sm2pull = pd.DataFrame(data=smooth(second_pull), columns=['strain', 'stress'])
        # sm2retract = pd.DataFrame(data=smooth(second_retract), columns=['strain', 'stress'])

        strain_num = []
        strain_num2 = []

        for i in range(len(Slist)):
            strain_num.append(strain(sm1pull, Slist[i]))
            strain_num2.append(strain(sm2pull, Slist[i]))
        S_group[g].loc[len(S_group[g])] = strain_num
        S_group2[g].loc[len(S_group2[g])] = strain_num2

        pla_defo = []

        pla_defo.append((ten_p(second_pull) - ori_p) / (5000 + ori_p) * 100)
        # pla_forplot = [ -x for x in pla_defo]
        p_defo.loc[len(p_defo)] = pla_defo
    #new_data = pd.concat([first_pull, first_retract, second_pull, second_retract], axis=1)

        #summary = pd.concat([summary, new_data],axis = 1)
        # ax = plt.subplot(111)
        # ax.plot(first_pull.strain *100, first_pull.stress, color = color[g], alpha = a)
        # ax.plot(first_retract.strain *100, first_retract.stress, color = color[g], linestyle = '--', alpha = a)
        # ax.plot(second_pull.strain *100, second_pull.stress, color = 'green')
        # ax.plot(second_retract.strain *100, second_retract.stress, color = 'green', linestyle = '--')
        # ax.plot(first_pull.stress, derv(sm1pull, inter), color = color[g])
        # ax.plot(first_retract.stress, derv(sm1retract, inter), color = color[g], linestyle = '--')

        # ax.set(ylabel='Stress (MPa)', xlabel = 'Strain (%)',title= Title + 'First cycle')
        # ax.set(ylabel='Load (gram)', xlabel='Strain (%)', title=Title + 'First cycle')
        # ax.set(ylabel='Modulus(MPa)', xlabel = 'Stress(MPa)',title= Title + '_first loading')
        # ax.grid(alpha=0.4, linestyle='--')

        # ax2 = plt.subplot(425)
        # ax2.plot(sm1pull.stress, derv(sm1pull, inter), color = color[g], alpha = a)
        # ax2.set(ylabel='Modulus (MPa)', xlabel = 'Stress (MPa)')
        # ax2.set(ylabel='Derivative (gram)', xlabel='Load (gram)')
        # ax2.grid(alpha=0.4, linestyle='--')

        # ax3 = plt.subplot(427)
        # ax3.plot(sm1retract.stress, derv(sm1retract, inter), color = color[g], linestyle = '--', alpha = a)
        # ax3.set(ylabel='Modulus (MPa)', xlabel = 'Stress (MPa)')
        # ax3.set(ylabel='Derivative (gram)', xlabel='Load (gram)')
        # ax3.grid(alpha=0.4, linestyle='--')


        # bx = plt.subplot(222)
        # bx.plot(second_pull.strain *100, second_pull.stress, color = color[g], alpha = a)
        # bx.plot(second_retract.strain *100, second_retract.stress, color = color[g], linestyle = '--', alpha = a)
        # bx.plot(first_pull.strain * 100, derv(sm1pull, inter), color = 'blue')
        # bx.plot(first_retract.strain *100, derv(sm1retract, inter), color = 'blue', linestyle = '--')
        # bx.plot(second_pull.stress, derv(sm2pull, inter), color = color[g])
        # bx.plot(second_retract.stress, derv(sm2retract, inter), color = color[g], linestyle = '--')
        # bx.set(ylabel='Stress(MPa)', xlabel='Strain(%)', title= Title + '_second cycle')
        # bx.set(ylabel='Stress (MPa)', xlabel='Strain (%)', title= Title + 'Second cycle')
        # bx.set(ylabel='Load (gram)', xlabel='Strain (%)', title=Title + 'Second cycle')
        # bx.grid(alpha=0.4, linestyle='--')

        # bx2 = plt.subplot(426)
        # bx2.plot(sm2pull.stress, derv(sm2pull, inter), color = color[g], alpha = a)
        # bx2.set(ylabel='Modulus (MPa)', xlabel = 'Stress (MPa)')
        # bx2.set(ylabel='Derivative (gram)', xlabel='Load (gram)')
        # bx2.grid(alpha=0.4, linestyle='--')

        # bx3 = plt.subplot(428)
        # bx3.plot(sm2retract.stress, derv(sm2retract, inter), color = color[g], linestyle = '--', alpha = a)
        # bx3.set(ylabel='Modulus (MPa)', xlabel = 'Stress (MPa)')
        # bx3.set(ylabel='Derivative (gram)', xlabel='Load (gram)')
        # bx3.grid(alpha=0.4, linestyle='--')


    # cx = plt.subplot(313)
        # cx.plot(first_pull.stress, derv(sm1pull, inter), color = 'blue')
        # cx.plot(first_retract.stress, derv(sm1retract, inter), color = 'blue', linestyle = '--')
        # cx.plot(second_pull.stress, derv(sm2pull, inter), color = 'green')
        # cx.plot(second_retract.stress, derv(sm2retract, inter), color = 'green', linestyle = '--')
        # cx.set(xlabel='Stress(MPa)', ylabel='Derivative(MPa)')
        # cx.grid(alpha=0.4, linestyle='--')

        # if n >= 3:
        #     break
        # plt.show()

        # bx = plt.subplot(121)
        # bx.plot(first_pull.strain, first_pull.stress, color = color[g], linestyle = style[g])
        # bx.set(ylabel='Stress (MPa)', xlabel='Strain')

        n += 1
        print(n)

    # bx = plt.subplot(121)
    # p_defo.to_csv(output[g])
    # P = np.mean(p_defo.defo)
    # print(p_defo, P)
    # P_std = np.std(p_defo.defo)
    # bx.scatter(1, P, color = color[g])


    g += 1

ax = plt.subplot(111)

y = Slist
ssF = []
ssN = []
sseF = []
sseN = []
s10 = []
s1 = []
s01 = []
se10 = []
se1 = []
se01 = []
#second pull
s10s = []
s1s = []
s01s = []
se10s = []
se1s = []
se01s = []

s10_4c = []
s1_4c = []
s01_4c = []
se10_4c = []
se1_4c = []
se01_4c = []
for i in range(len(Slist)):
    # ssF.append(np.mean(ssFresh.iloc[:, i]))
    # ssN.append(np.mean(ssNext.iloc[:, i]))
    # sseF.append(np.std(ssFresh.iloc[:, i])/(len(ssFresh[0])) ** 1/2)
    # sseN.append(np.std(ssNext.iloc[:, i])/(len(ssNext[0])) ** 1/2)
    s10.append(np.mean(R10c.iloc[:, i]) *100)
    s1.append(np.mean(R1c.iloc[:, i]) *100)
    s01.append(np.mean(R01c.iloc[:, i]) *100)
    s10s.append(np.mean(R10c2.iloc[:, i]) *100)
    s1s.append(np.mean(R1c2.iloc[:, i]) *100)
    s01s.append(np.mean(R01c2.iloc[:, i]) *100)

    # se10.append(np.std(R10c.iloc[:, i]))
    # se1.append(np.std(R1c.iloc[:, i]))
    # se01.append(np.std(R01c.iloc[:, i]))

    # s10_4c.append(np.mean(R10c_4c.iloc[:, i]))
    # s1_4c.append(np.mean(R1c_4c.iloc[:, i]))
    # s01_4c.append(np.mean(R01c_4c.iloc[:, i]))
    # se10_4c.append(np.std(R10c_4c.iloc[:, i]))
    # se1_4c.append(np.std(R1c_4c.iloc[:, i]))
    # se01_4c.append(np.std(R01c_4c.iloc[:, i]))

    se10.append((np.std(R10c.iloc[:, i]) *100)/(len(R10c[0]) ** 1/2))
    se1.append((np.std(R1c.iloc[:, i]) *100)/(len(R1c[0]) ** 1/2))
    se01.append((np.std(R01c.iloc[:, i]) *100)/(len(R01c[0]) ** 1/2))
    se10s.append((np.std(R10c2.iloc[:, i]) *100)/(len(R10c2[0]) ** 1/2))
    se1s.append((np.std(R1c2.iloc[:, i]) *100)/(len(R1c2[0]) ** 1/2))
    se01s.append((np.std(R01c2.iloc[:, i]) *100)/(len(R01c2[0]) ** 1/2))

    # se10.append((np.std(R10c.iloc[:, i]) *100))
    # se1.append((np.std(R1c.iloc[:, i]) *100))
    # se01.append((np.std(R01c.iloc[:, i]) *100))
    # se10s.append((np.std(R10c2.iloc[:, i]) *100))
    # se1s.append((np.std(R1c2.iloc[:, i]) *100))
    # se01s.append((np.std(R01c2.iloc[:, i]) *100))

    #
    # se10_4c.append(np.std(R10c_4c.iloc[:, i])/(len(R10c_4c[0])) ** 1/2)
    # se1_4c.append(np.std(R1c_4c.iloc[:, i])/(len(R1c_4c[0])) ** 1/2)
    # se01_4c.append(np.std(R01c_4c.iloc[:, i])/(len(R01c_4c[0])) ** 1/2)
# ax.errorbar(ssF , y, xerr = sseF, color= 'blue', alpha = 0.6)
# ax.errorbar(ssN , y, xerr = sseN, color='red', alpha = 0.6)
# ax.errorbar(s10 , y, xerr = se10, color= 'orange', alpha = 0.6)
# ax.errorbar(s1 , y, xerr = se1, color='red', alpha = 0.6)
# ax.errorbar(s01 , y, xerr = se01, color='purple', alpha = 0.6)
ax.errorbar(s10s , y, xerr = se10s, color= 'orange', alpha = 0.6)
ax.errorbar(s1s , y, xerr = se1s, color='red', alpha = 0.6)
ax.errorbar(s01s , y, xerr = se01s, color='purple', alpha = 0.6)

#
# ax.errorbar(s10_4c , y, xerr = se10_4c, color= 'orange', alpha = 0.6, linestyle = 'dashed')
# ax.errorbar(s1_4c , y, xerr = se1_4c, color='red', alpha = 0.6, linestyle = 'dashed')
# ax.errorbar(s01_4c , y, xerr = se01_4c, color='purple', alpha = 0.6, linestyle = 'dashed')

# ax.errorbar(s10 , y, color= 'orange', alpha = 0.9)
# ax.errorbar(s1 , y, color='blue', alpha = 0.9)
# ax.errorbar(s01 , y, color='green', alpha = 0.9)
#
# ax.errorbar(s10_4c , y, color= 'red', alpha = 0.9)
# ax.errorbar(s1_4c , y, color='black', alpha = 0.9)
# ax.errorbar(s01_4c , y, color='purple', alpha = 0.9)

ax.set(ylabel='Stress (MPa)', xlabel = 'Strain', ylim = [-0.5, 5.5])


# orange_line = mlines.Line2D([], [], color = 'orange', label = '10 mm/min')
# blue_line = mlines.Line2D([], [], color = 'blue', label = '1 mm/min')
# green_line = mlines.Line2D([], [], color = 'green', label = '0.1 mm/min')
# red_line = mlines.Line2D([], [], color = 'red', label = '10 mm/min_4c')
# black_line = mlines.Line2D([], [], color = 'black', label = '1 mm/min_4c')
# purple_line = mlines.Line2D([], [], color = 'purple', label = '0.1 mm/min_4c')
# orange_line = mlines.Line2D([], [], color = 'blue', label = 'Fresh')
# blue_line = mlines.Line2D([], [], color = 'red', label = 'Next day')
orange_line = mlines.Line2D([], [], color = 'orange', label = '10 mm/min')
blue_line = mlines.Line2D([], [], color = 'red', label = '1 mm/min')
green_line = mlines.Line2D([], [], color = 'purple', label = '0.1 mm/min')
# red_line = mlines.Line2D([], [], color = 'orange', label = '10 mm/min_4c', linestyle = 'dashed')
# black_line = mlines.Line2D([], [], color = 'red', label = '1 mm/min_4c', linestyle = 'dashed')
# purple_line = mlines.Line2D([], [], color = 'purple', label = '0.1 mm/min_4c', linestyle = 'dashed')
ax.legend(handles = [orange_line, blue_line, green_line], loc = 2, fontsize=10, frameon=False)

# blue_line = mlines.Line2D([], [], color = 'C1', label = 'Retract compliance', linestyle = '--')
# ax.legend(handles = [orange_line, blue_line, green_line, red_line, black_line, purple_line], loc = 2, fontsize = 8)
# bx.legend(handles = [orange_line, blue_line, green_line, red_line, black_line, purple_line], loc = 2, fontsize = 8)

# ax.legend(handles = [black_line, blue_line, green_line, orange_line], loc = 2, fontsize = 8)
# bx.legend(handles = [black_line, blue_line, green_line, orange_line], loc = 2, fontsize = 8)

plt.gcf().subplots_adjust(bottom=0.15, right= 0.83)
plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_4)_secondpull.pdf', transparent=True)
plt.show()

