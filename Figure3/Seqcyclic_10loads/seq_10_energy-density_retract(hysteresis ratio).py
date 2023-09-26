# script: Study progression of elasticity and plasticity with strain energy density method
# output value of T & E & P strain energy density and their incremental strain energy density;
# plot: line-plot for exact value, stack bar chart for incremental value

# import necessary package
import pandas as pd
import numpy as np
import os
import glob
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

# User input
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'

group = [folder2, Onion3, Onion5, Onion6]

thickness = 7
width = 3
smf = 0.04

# for load sequence of 4-8-12-16-20
loadx = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
# for load sequence of 8-12-16-20-24
# loadx = [8, 12, 16, 20, 24]

# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

# define a function that find the index where peels are under tension ( constant > 0.1g )
def ten_ind(pull, load):
    for i in range(len(pull) - 5):
        if (pull.load[i] > load) & (pull.load[i + 1] > load) & (pull.load[i + 2] > load) & \
                (pull.load[i + 3] > load) & (pull.load[i + 4] > load):
            index = i
            return index
        elif i > len(pull) - 7:
            print('no fit')

# define a function that find the position where peels are under tension ( constant > 0.1g )
def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

# define a function that calculate the stress and strain of the curve
def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / 5000

def smooth(pull, smf):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= smf)
    return z

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

# define a function to calculate the area for total pull
def area_t(pull):
    area = 0
    for i in range(len(pull) - 1):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    return area

# calculate the area to certain load
def area_f(pull, load):
    area = 0
    for i in range(ten_ind(pull, load)):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    return area

# calculate the area according to the previous pull max strain
def area_s(pull, pre_pull):
    area = 0
    cutted_data = pull[pull.strain < np.max(pre_pull.strain)].reset_index(drop=True)
    for i in range(len(cutted_data)):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    return area

# create data frame for results from each property
t_energ = pd.DataFrame(columns=loadx)
e_energ = pd.DataFrame(columns=loadx)
p_energ = pd.DataFrame(columns=loadx)
store_energ = pd.DataFrame(columns=loadx)
hyst_energy = pd.DataFrame(columns=loadx)

ela_ratio_DF = pd.DataFrame(columns=loadx)
pla_ratio_DF = pd.DataFrame(columns=loadx)
sto_ratio_DF = pd.DataFrame(columns=loadx)

t_inc_energ = pd.DataFrame(columns=loadx)
e_inc_energ = pd.DataFrame(columns=loadx)
p_inc_energ = pd.DataFrame(columns=loadx)
h_inc_energ = pd.DataFrame(columns=loadx)
s_inc_energ = pd.DataFrame(columns=loadx)

n = 0
# extract data from each pull
# Use if multiple files need analysis
g = 0
h = ['use default', '/']
for folder in group:

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
            n += 1

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull, eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

        ### normalize curve
        for i in range(len(pull)):
            norm(pull[i])
            norm(retract[i])

        # target for seuqntial Instron
        target = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20]
        # loop for seuqntial Instron
        for i in range(len(pull)):
            retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace = True)

            retract[i].reset_index(inplace=True)
            rmindex = []
            purify_revrs(retract[i])

            retract[i].drop(retract[i].loc[rmindex].index, inplace=True)
            retract[i] = retract[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)

        sm1pull = pd.DataFrame(data=smooth(first_pull, smf), columns=['strain', 'stress'])
        sm2pull = pd.DataFrame(data=smooth(second_pull, smf), columns=['strain', 'stress'])
        sm3pull = pd.DataFrame(data=smooth(third_pull, smf), columns=['strain', 'stress'])
        sm4pull = pd.DataFrame(data=smooth(fourth_pull, smf), columns=['strain', 'stress'])
        sm5pull = pd.DataFrame(data=smooth(fifth_pull, smf), columns=['strain', 'stress'])
        sm7pull = pd.DataFrame(data=smooth(seventh_pull, smf), columns=['strain', 'stress'])
        sm8pull = pd.DataFrame(data=smooth(eighth_pull, smf), columns=['strain', 'stress'])
        sm9pull = pd.DataFrame(data=smooth(ninth_pull, smf), columns=['strain', 'stress'])
        sm10pull = pd.DataFrame(data=smooth(tenth_pull, smf), columns=['strain', 'stress'])
        sm11pull = pd.DataFrame(data=smooth(eleventh_pull, smf), columns=['strain', 'stress'])
        sm12pull = pd.DataFrame(data=smooth(twlveth_pull, smf), columns=['strain', 'stress'])


        sm1retract = pd.DataFrame(data=smooth(first_retract, smf), columns=['strain', 'stress'])
        sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
        sm3retract = pd.DataFrame(data=smooth(third_retract, smf), columns=['strain', 'stress'])
        sm4retract = pd.DataFrame(data=smooth(fourth_retract, smf), columns=['strain', 'stress'])
        sm5retract = pd.DataFrame(data=smooth(fifth_retract, smf), columns=['strain', 'stress'])
        sm7retract = pd.DataFrame(data=smooth(seventh_retract, smf), columns=['strain', 'stress'])
        sm8retract = pd.DataFrame(data=smooth(eighth_retract, smf), columns=['strain', 'stress'])
        sm9retract = pd.DataFrame(data=smooth(ninth_retract, smf), columns=['strain', 'stress'])
        sm10retract = pd.DataFrame(data=smooth(tenth_retract, smf), columns=['strain', 'stress'])
        sm11retract = pd.DataFrame(data=smooth(eleventh_retract, smf), columns=['strain', 'stress'])
        sm12retract = pd.DataFrame(data=smooth(twlveth_retract, smf), columns=['strain', 'stress'])

        smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm7pull, sm8pull, sm9pull, sm10pull, sm11pull, sm12pull]
        smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm7retract, sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]

        for i in range(len(smpull)):
            smpull[i]['load'] = pull[i].load
            smretract[i]['load'] = retract[i].load

        loadseq = [4, 8, 12, 16, 20]

        # incremental plastic strain energy density
        inc_p = []
        for i in range(len(smpull)-1):
            inc_p.append(area_t(smpull[i]) - area_s(smpull[i+1], smpull[i]))

        p_inc_energ.loc[len(p_inc_energ)] = inc_p

        # for sequential Instron - plastic strain energy density
        pla_energ = []
        ener = 0
        for i in range(len(inc_p)):
            ener += inc_p[i]
            pla_energ.append(ener)
        p_energ.loc[len(p_energ)] = pla_energ

        # date reporter
        print(p_energ)

        # for sequential Instron - elastic strain energy density
        ela_energ = []
        for i in range(len(smpull)-1):
            ela_energ.append(area_s(smpull[i+1], smpull[i]))
        e_energ.loc[len(e_energ)] = ela_energ

        # for sequential Instron - total strain energy density
        total_energ = []
        for i in range(len(ela_energ)):
            total_energ.append(ela_energ[i] + pla_energ[i])

        # stored energy
        retra = []
        for i in range(len(smretract)-1):
            retra.append(area_t(smretract[i]))
        store_energ.loc[len(store_energ)] = retra

        # hysteresis dissipated energy
        hys = []
        for i in range(len(loadx)):
            hys.append(ela_energ[i] - retra[i])
        hyst_energy.loc[len(hyst_energy)] = hys

    # add data of this file to the results data frame
        t_energ.loc[len(t_energ)] = total_energ

    # incremental total strain energy density
        inc_t = [total_energ[0]]
        for i in range(len(ela_energ)):
            if i > 0:
                inc_t.append(total_energ[i] - total_energ[i - 1])
        t_inc_energ.loc[len(t_inc_energ)] = inc_t

    ## ratio calculatino
    # plastic energy dissipation
        pla_ratio = []
        for i in range(len(loadx)):
            pla_ratio.append(pla_energ[i]/total_energ[i] *100)
        pla_ratio_DF.loc[len(pla_ratio_DF)] = pla_ratio

        ela_ratio = []
        for i in range(len(loadx)):
            ela_ratio.append(hys[i]/total_energ[i] *100)
        ela_ratio_DF.loc[len(ela_ratio_DF)] = ela_ratio

        sto_ratio = []
        for i in range(len(loadx)):
            sto_ratio.append(retra[i]/total_energ[i] *100)
        sto_ratio_DF.loc[len(sto_ratio_DF)] = sto_ratio

# plot as stress
x = [a * 0.098 /3/7/0.01 for a in loadx]

C_p = []
C_h = []
C_s = []
Ce_p = []
Ce_h = []
Ce_s = []
for i in range(len(loadx)):
    C_p.append(np.mean(pla_ratio_DF.iloc[:, i]))
    C_h.append(np.mean(ela_ratio_DF.iloc[:, i]))
    C_s.append(np.mean(sto_ratio_DF.iloc[:, i]))
    Ce_p.append(np.std(pla_ratio_DF.iloc[:, i]))
    Ce_h.append(np.std(ela_ratio_DF.iloc[:, i]))
    Ce_s.append(np.std(sto_ratio_DF.iloc[:, i]))

# bar-plot with error bar
bx = plt.subplot(111)
xlabel = [0, 1, 2, 3, 4, 5, 6 , 7, 8, 9, 10]
lw = 2
ms = 4
bx.errorbar(x, C_p, yerr = Ce_p, marker = '8',color='C9', alpha = 0.6, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
bx.errorbar(x, C_h, yerr = Ce_h,  marker = '8',color='C3', alpha = 0.6, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
bx.errorbar(x, C_s, yerr = Ce_s,  marker = '8',color='C2', alpha = 0.6, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
blue_line = mlines.Line2D([], [], color = 'C2', label = 'Stored energy' )
green_line = mlines.Line2D([], [], color = 'C3', label = 'Dissipated energy (elastic hysteresis)')
orange_line = mlines.Line2D([], [], color = 'C9', label = 'Dissipated energy(plasticity)')
bx.legend(handles = [blue_line, green_line, orange_line], loc = 2, fontsize=10, frameon = False)

bx.set_xticks(xlabel)
bx.set_xticklabels(xlabel)
bx.set_yticks(range(0, 80, 10))
bx.set(ylim = [-10, 80], xlabel='Stress (MPa)', ylabel='Percentage (%)')

plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
# plt.savefig(''f'{output}/Energy_density_ratio.pdf', transparent = True)
plt.show()

