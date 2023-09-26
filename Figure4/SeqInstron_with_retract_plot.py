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

thickness = 7
width = 3
smf = 0.04
intv = 5
l = -0.5
r = 45

folder = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/Plot'



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

def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return (len(pull) -1)

# define a function that calculate the stress and strain of the curve
def norm(pull):
    pull['force_N'] = pull.load * 0.0098
    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)


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

def ten_ind(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            return i

color = ['green', 'orange', 'blue', 'black']
# color = ['orange', 'blue', 'black']

n = 0
for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
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


    pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
    retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]

    # target for seuqntial Instron
    target = [4, 8, 12, 16, 20, 20]
    ori_p = ten_p(first_pull)
    # loop for seuqntial Instron
    for i in range(len(pull)):
        # retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)
        retract[i].drop(retract[i].loc[0:ten_ind(retract[i])].index, inplace=True)
        pull[i].drop(pull[i].loc[0:ten_ind(pull[i])].index, inplace = True)

        retract[i].reset_index(inplace=True)
        rmindex = []
        purify_revrs(retract[i])
        retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

        retract[i].reset_index(inplace=True)
        pull[i].reset_index(inplace=True)

    ### normalize curve
    for i in range(len(pull)):
        norm(pull[i])
        norm(retract[i])

    sm1pull = pd.DataFrame(data = smooth(first_pull,smf),columns=['strain','stress'])
    sm2pull = pd.DataFrame(data = smooth(second_pull,smf),columns=['strain','stress'])
    sm3pull = pd.DataFrame(data = smooth(third_pull,smf),columns=['strain','stress'])
    sm4pull = pd.DataFrame(data = smooth(fourth_pull,smf),columns=['strain','stress'])
    sm5pull = pd.DataFrame(data = smooth(fifth_pull,smf),columns=['strain','stress'])
    sm6pull = pd.DataFrame(data = smooth(sixth_pull,smf),columns=['strain','stress'])

    sm1retract = pd.DataFrame(data=smooth(first_retract,smf), columns=['strain','stress'])
    sm2retract = pd.DataFrame(data=smooth(second_retract,smf), columns=['strain','stress'])
    sm3retract = pd.DataFrame(data = smooth(third_retract,smf),columns=['strain','stress'])
    sm4retract = pd.DataFrame(data = smooth(fourth_retract,smf),columns=['strain','stress'])
    sm5retract = pd.DataFrame(data = smooth(fifth_retract,smf),columns=['strain','stress'])
    sm6retract = pd.DataFrame(data = smooth(sixth_retract,smf),columns=['strain','stress'])
    smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull]
    smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract]

    first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(
        drop=True)
    first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[
                    ::-1].reset_index(drop=True)
    second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(
        drop=True)
    second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[
                     ::-1].reset_index(drop=True)
    third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(
        drop=True)
    third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[
                    ::-1].reset_index(drop=True)
    fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(
        drop=True)
    fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[
                     ::-1].reset_index(drop=True)
    fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(
        drop=True)
    fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[
                    ::-1].reset_index(drop=True)
    sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(
        drop=True)
    sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)

    pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
    retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]

    for i in range(len(pull)):
        norm(pull[i])
        norm(retract[i])

    map = 'Paired'
    ax = plt.subplot(111)
    # for i in range(len(pull)):
# for only first two cycle
    for i in range(len(smretract)):
        # unsmoothed data
        ax.plot(pull[i].strain *100, pull[i].stress, color = color[n])
        ax.plot(retract[i].strain *100, retract[i].stress, color = color[n], linestyle = '--')

    ax.set(ylabel = 'Stress (MPa)', xlabel = 'Strain (%)')

    n += 1
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(12.8, 5.5)
green_line = mlines.Line2D([], [], color='green', label='Control')
orange_line = mlines.Line2D([], [], color='orange', label='40% PEG8k')
blue_line = mlines.Line2D([], [], color='blue', label='40% PEG8k rehydrated')
ax.legend(handles=[green_line, orange_line, blue_line], loc=2, fontsize=10, frameon=False)

# plt.savefig(''f'{output}/Con_PEG_rehydra_seq_plot.pdf', transparent = True)

plt.show()





