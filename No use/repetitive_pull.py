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
fig.set_size_inches(12.8, 4.8)

thickness = 7
width = 3
smf = 0.04
intv = 5
l = -0.5
r = 45


repetitive10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/repetitive pulling_8.1.20/10g'
repetitive10_2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/repetitive_pulling_10g_3mpm_10.6.20'

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


# path = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Nitro_glove_dif_rate_20g/Instron_with_Retract__03_18_2020__15_06_56_SHORT.csv'

#for file in sorted(glob.glob(os.path.join(path, '*SHORT.csv'))):

color = ['blue','green', 'orange', 'red', 'purple', 'black']
n = 1
for file in sorted(glob.glob(os.path.join(repetitive10_2, '*SHORT.csv'))):
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
    ax = plt.subplot(121)
    # for i in range(len(pull)):
# for only first two cycle
    for i in range(len(smretract)):
        # unsmoothed data
        ax.plot(pull[i].strain *100, pull[i].stress, color = color[i])
        ax.plot(retract[i].strain *100, retract[i].stress, color = color[i], linestyle = '--')

    ax.set(ylabel = 'Stress(MPa)', xlabel = 'Strain(%)')

    # setting for gram y axis
    ax0y = ax.twinx()
    ax0y.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
    B = -1 * 9.8
    T = 12 * 9.8
    ax0y.set(ylim=[B, T])
    SB = B * 0.001 / (thickness * width * 0.001)
    ST = T * 0.001 / (thickness * width * 0.001)
    ax.set(ylim=[SB, ST])

    #
    bx = plt.subplot(122)
    #
    for j in range(len(smpull)):
        bx.plot(smpull[j].strain *100, derv(smpull[j], intv) , color = color[j])

    bx.set(ylabel = 'Modulus (MPa)', xlabel = 'Strain(%)')
    # bx.set(ylim = [10**-2, 10**1.05],ylabel = 'Stress(MPa)', xlabel = 'Strain(%)')
    # bx.grid(alpha=0.4, linestyle='--')
    # bx.set_yscale('log')

    # dx = plt.subplot(313)
    for j in range(len(smretract)):
        bx.plot(smretract[j].strain *100, derv(smretract[j], intv) , color = color[j], linestyle = '--')

    #     dx.scatter(smpull[j].strain * 100, dervuni(np.log(smpull[j].stress), smpull[j].strain, intv),
    #                c=derv(smpull[j], intv), cmap=map, s=2)
    #     dx.scatter(smretract[j].strain * 100, dervuni(np.log(smretract[j].stress), smretract[j].strain, intv),
    #                c=derv(smretract[j], intv), cmap=map, s=2)
        # dx.plot(smpull[j].strain * 100, dervuni(np.log(smpull[j].stress), smpull[j].strain, intv), color = color[j])
        # dx.plot(smretract[j].strain * 100, dervuni(np.log(smretract[j].stress), smretract[j].strain, intv), color = color[j],linestyle = '--')
        # dx.plot(smpull[j].strain *100, smpull[j].stress, color = color[j])
        # dx.scatter(smpull[j].strain * 100, smpull[j].stress, c=derv(smpull[j], intv), cmap=map, s=1)
        # dx.scatter(smretract[j].strain * 100, smretract[j].stress, c=derv(smretract[j], intv), cmap=map, s=1)
    #
    # dx.set(ylabel = 'Derivative (MPa)', xlabel = 'Strain(%)')
    # dx.grid(alpha=0.4, linestyle='--')
    # dx.set(ylim = [10**-2, 10**1.05], ylabel = 'Stress (MPa)', xlabel = 'Strain(%)')
    # dx.set_yscale('log')

    # dx.set(ylim = [-10,300])
    # bx2 = bx.twiny()
    # bx2.set(xlim = [l*0.01*5000/400*1.18, r*0.01*5000/400*1.18], xlabel = 'recording time(s)')
    #
    # cx = plt.subplot(312)
    #
    # for j in range(len(smpull)):
    #     cx.scatter(smpull[j].strain *100, derv(smpull[j], intv), c=derv(smpull[j], intv), cmap = map, s =2)
        # cx.scatter(smretract[j].strain *100, derv(smretract[j], intv), c=derv(smretract[j], intv), cmap = map, s =2)
        # cx.plot(smpull[j].strain *100, derv(smpull[j], intv), color = color[j])
        # cx.plot(smretract[j].strain *100, derv(smretract[j], intv), color = color[j], linestyle = '--')
        # cx.scatter(smretract[j].stress , derv(smretract[j], intv) , c=derv(smretract[j], intv), cmap = map, s =2)

    # cx.set(ylim = [-5,130], xlabel = 'Strain(%)', ylabel = 'Derivative(MPa)')
    # cx.grid(alpha = 0.4, linestyle = '--')



# when we plot loading and unloading curve separatly
#     dx = plt.subplot(313)




    # for j in range(len(smpull)):
        # dx.plot(smpull[j].strain * 100, derv(smpull[j], intv), color=color[j])
        # dx.plot(smretract[j].strain * 100, derv(smretract[j], intv), color=color[j], linestyle='--')
        # bx.scatter(smretract[j].strain , derv(smretract[j], intv) , c=derv(smretract[j], intv), cmap = map, s =2)

    # dx.set(xlim = [l,r], ylabel='Derivative(MPa)', xlabel='Strain(%)')
    # dx.grid(alpha=0.4, linestyle='--')
    #
    # dx2 = dx.twiny()
    # dx2.set_xlim([l*0.01*5000/400*1.18, r*0.01*5000/400*1.18])
    # dx2.set_xlabel('recording time(s)')

    # ax = plt.subplot(111)
    blue_line = mlines.Line2D([], [], color = 'blue', label = '1˚cycle' )
    green_line = mlines.Line2D([], [], color = 'green', label = '2˚cycle' )
    orange_line = mlines.Line2D([], [], color = 'orange', label = '3˚cycle' )
    red_line = mlines.Line2D([], [], color = 'red', label = '4˚cycle' )
    purple_line = mlines.Line2D([], [], color = 'purple', label = '5˚cycle' )
    black_line = mlines.Line2D([], [], color = 'black', label = '6˚cycle' )
    # ax.set(xlabel = 'strain', ylabel = 'Force (N)', title = 'Nitro-glove cyclic pulling_10mm/min')
    ax.legend(handles = [blue_line, green_line, orange_line, red_line, purple_line, black_line], loc = 2, fontsize=10, frameon = False)
    bx.legend(handles = [blue_line, green_line, orange_line, red_line, purple_line, black_line], loc = 2, fontsize=10, frameon = False)

    n += 1
    plt.gcf().subplots_adjust(bottom=0.2, wspace= 0.4)
    fig = plt.gcf()
    size = fig.get_size_inches()
    print(size)
    # plt.savefig(''f'{output}/Repetitive_pulling.pdf', transparent = True)

    plt.show()

# ax = plt.subplot(111)
# blue_line = mlines.Line2D([], [], color = 'blue', label = '1˚cycle' )
# green_line = mlines.Line2D([], [], color = 'green', label = '2˚cycle' )
# orange_line = mlines.Line2D([], [], color = 'orange', label = '3˚cycle' )
# black_line = mlines.Line2D([], [], color = 'black', label = '4˚cycle' )
# ax.set(xlabel = 'strain', ylabel = 'Force (N)', title = 'Nitro-glove cyclic pulling_10mm/min')
# ax.legend(handles = [blue_line, green_line, orange_line, black_line])




