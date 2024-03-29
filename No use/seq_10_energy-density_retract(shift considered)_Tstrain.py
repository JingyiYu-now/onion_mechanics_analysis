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
from numpy import trapz
import statsmodels.api as sm
import matplotlib as mpl
mpl.use('macosx')
# User input

# PEG Ca2+ sample
# Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/buffer_con'
# Con_Ca = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/Ca2+_con'
# PEG_con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/buffer+PEG'
# PEG_Ca = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/Ca2+_PEG'

#Onion experiment with rehydration
Con1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/buffer'
PEG1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/PEG_2h'
Reh1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/PEG re-hydrated'

Con2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/Control'
PEG2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/PEG2h'
Reh2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/PEG_2h_rehydrated'

Con3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/Control'
PEG3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/PEG_2h'
Reh3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/PEG_2h_rehydration'


output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/PEG'
# group = [ folder2, Onion3, Onion5]
# group = [Con4, PEG4]
# group = [ folder2]
Ucon = [Con1, Con2, Con3]
UPEG = [PEG1, PEG2, PEG3]
UReH = [Reh1, Reh2, Reh3]
group = [Ucon, UPEG]

name = 'Onion5'
# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output1 = ''f'{folder}/summary1.csv'
# output2 = ''f'{folder}/summary_inc_1.csv'
# output = ''f'{folder3}/total_chart_1.csv'

thickness = 7
width = 3
smf = 0.04
marker = ['8', 's']
# for load sequence of 4-8-12-16-20
strainT = [5, 10, 15, 20, 25]
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

# Define a function that find the position where peels are under tension ( constant > 0.1g )
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
    # strain calculated based on length at the beginning of each pull
    # ori_p = ten_p(pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / 5000
    # true stress & strain
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p + 4500) )
    # pull['stress'] = pull.force_N/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)


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


# calculate the area according to the previous pull max strain
def area_s(pull, pre_pull):
    area = 0
    cutted_data = pull[pull.strain < np.max(pre_pull.strain)].reset_index(drop=True)
    for i in range(len(cutted_data)):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    return area

# create data frame for results from each property


n = 0
# extract data from each pull
# Use if multiple files need analysis
g = 0
h = ['use default', '/']
ls = ['-', '--']
for unit in group:
    t_energ = pd.DataFrame(columns=strainT)
    e_energ = pd.DataFrame(columns=strainT)
    p_energ = pd.DataFrame(columns=strainT)
    store_energ = pd.DataFrame(columns=strainT)
    hyst_energy = pd.DataFrame(columns=strainT)
    t_inc_energ = pd.DataFrame(columns=strainT)
    e_inc_energ = pd.DataFrame(columns=strainT)
    p_inc_energ = pd.DataFrame(columns=strainT)
    h_inc_energ = pd.DataFrame(columns=strainT)
    s_inc_energ = pd.DataFrame(columns=strainT)

    for folder in unit:

        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        # if n%2 == 0:
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
            # continue

            # if n%2 == 1:
            #     print(n)
            #     whole_df = pd.read_table(file,
            #                          delimiter=',',
            #                          header=None,
            #                          names=range(0, 2),
            #                          )
            #     whole_df.columns = ['position', 'load']
            #
            #
            #     seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
            #     seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            #     eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
            #     eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            #     ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(drop=True)
            #     ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            #     tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(drop=True)
            #     tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            #     eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(drop=True)
            #     eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            #     twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(drop=True)
            #     twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)
            #     n += 1
            #
            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]


            ### normalize curve

            for i in range(len(pull)):
                norm(pull[i])
                norm(retract[i])

            # target for seuqntial Instron
            # target = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20]
            # loop for seuqntial Instron
            for i in range(len(pull)):
                # retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace = True)
                # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
                # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

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
            sm6pull = pd.DataFrame(data=smooth(sixth_pull, smf), columns=['strain', 'stress'])

        # sm7pull = pd.DataFrame(data=smooth(seventh_pull, smf), columns=['strain', 'stress'])
            # sm8pull = pd.DataFrame(data=smooth(eighth_pull, smf), columns=['strain', 'stress'])
            # sm9pull = pd.DataFrame(data=smooth(ninth_pull, smf), columns=['strain', 'stress'])
            # sm10pull = pd.DataFrame(data=smooth(tenth_pull, smf), columns=['strain', 'stress'])
            # sm11pull = pd.DataFrame(data=smooth(eleventh_pull, smf), columns=['strain', 'stress'])
            # sm12pull = pd.DataFrame(data=smooth(twlveth_pull, smf), columns=['strain', 'stress'])


            sm1retract = pd.DataFrame(data=smooth(first_retract, smf), columns=['strain', 'stress'])
            sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
            sm3retract = pd.DataFrame(data=smooth(third_retract, smf), columns=['strain', 'stress'])
            sm4retract = pd.DataFrame(data=smooth(fourth_retract, smf), columns=['strain', 'stress'])
            sm5retract = pd.DataFrame(data=smooth(fifth_retract, smf), columns=['strain', 'stress'])
            sm6retract = pd.DataFrame(data=smooth(sixth_retract, smf), columns=['strain', 'stress'])

        # sm7retract = pd.DataFrame(data=smooth(seventh_retract, smf), columns=['strain', 'stress'])
            # sm8retract = pd.DataFrame(data=smooth(eighth_retract, smf), columns=['strain', 'stress'])
            # sm9retract = pd.DataFrame(data=smooth(ninth_retract, smf), columns=['strain', 'stress'])
            # sm10retract = pd.DataFrame(data=smooth(tenth_retract, smf), columns=['strain', 'stress'])
            # sm11retract = pd.DataFrame(data=smooth(eleventh_retract, smf), columns=['strain', 'stress'])
            # sm12retract = pd.DataFrame(data=smooth(twlveth_retract, smf), columns=['strain', 'stress'])

            smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull]
            smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract]



            for i in range(len(smpull)):
                smpull[i]['load'] = pull[i].load
                smretract[i]['load'] = retract[i].load



            # incremental plastic strain energy density
            inc_p = []
            for i in range(len(smpull)-2):
                inc_p.append(area_t(smpull[i]) - area_s(smpull[i+1], smpull[i]))
                print(inc_p)
            inc_p.append(area_t(smpull[len(smpull)-2]) - area_t(smpull[len(smpull)-1]))
            print(inc_p)
            # inc_p.append(area_t(sm1pull) - area_f(sm2pull, loadx[0]))
            # inc_p.append(area_t(sm2pull) - area_f(sm3pull, loadx[1]))
            # inc_p.append(area_t(sm3pull) - area_f(sm4pull, loadx[2]))
            # inc_p.append(area_t(sm4pull) - area_f(sm5pull, loadx[3]))
            # inc_p.append(area_t(sm5pull) - area_t(sm6pull))
            p_inc_energ.loc[len(p_inc_energ)] = inc_p

            # for sequential Instron - plastic strain energy density
            pla_energ = []
            ener = 0
            for i in range(len(inc_p)):
                ener += inc_p[i]
                pla_energ.append(ener)
            p_energ.loc[len(p_energ)] = pla_energ

            # for sequential Instron - elastic strain energy density
            ela_energ = []
            for i in range(len(smpull)-2):
                ela_energ.append(area_s(smpull[i+1], smpull[i]))
            ela_energ.append(area_t(smpull[len(smpull)-1]))
            e_energ.loc[len(e_energ)] = ela_energ
            # ela_energ.append(area_f(sm2pull, loadx[0]))
            # ela_energ.append(area_f(sm3pull, loadx[1]))
            # ela_energ.append(area_f(sm4pull, loadx[2]))
            # ela_energ.append(area_f(sm5pull, loadx[3]))
            # ela_energ.append(area_t(sm12pull))


            # incremental elastic strain energy density
            # inc_e = []
            # for i in range(len(smpull)-2):
            #     if i == 0:
            #         inc_e.append(area_f(smpull[i+1], loadx[0]))
            #     else:
            #         inc_e.append(area_f(smpull[i+2], loadx[1]) - area_f(smpull[i+1], loadx[0]))
            # inc_e.append(area_f(sm3pull, loadx[1]) - area_f(sm2pull, loadx[0]))
            # inc_e.append(area_f(sm4pull, loadx[2]) - area_f(sm3pull, loadx[1]))
            # inc_e.append(area_f(sm5pull, loadx[3]) - area_f(sm4pull, loadx[2]))
            # inc_e.append(area_t(sm6pull) - area_f(sm5pull, loadx[3]))
            # e_inc_energ.loc[len(e_inc_energ)] = inc_e

            # for sequential Instron - total strain energy density
            total_energ = []
            for i in range(len(ela_energ)):
                total_energ.append(ela_energ[i] + pla_energ[i])

            # stored energy
            retra = []
            for i in range(len(smretract)-1):
                retra.append(area_t(smretract[i]))
            store_energ.loc[len(store_energ)] = retra

            # incremental stored energy
            # sto_inc = []
            # for i in range(len(loadx)):
            #     if i == 0:
            #         sto_inc.append(retra[i])
            #     else:
            #         sto_inc.append(retra[i] - retra[i - 1])
            # s_inc_energ.loc[len(s_inc_energ)] = sto_inc

            # hysteresis dissipated energy
            hys = []
            for i in range(len(strainT)):
                hys.append(ela_energ[i] - retra[i])
            hyst_energy.loc[len(hyst_energy)] = hys

            # incremental hysteresis energy
            # hys_inc = []
            # for i in range(len(loadx)):
            #     if i == 0:
            #         hys_inc.append(hys[i])
            #     else:
            #         hys_inc.append(hys[i] - hys[i-1])
            # h_inc_energ.loc[len(h_inc_energ)] = hys_inc

        # add data of this file to the results data frame
            t_energ.loc[len(t_energ)] = total_energ

        # incremental total strain energy density
        #     inc_t = [total_energ[0]]
        #     for i in range(len(ela_energ)):
        #         if i > 0:
        #             inc_t.append(total_energ[i] - total_energ[i - 1])
        #     t_inc_energ.loc[len(t_inc_energ)] = inc_t






    ###### results plotting
    ##### set x
        # nb = 1
        # t = 2
        # d = 10
        # w = 0.8
        # xc = [t * x + w * nb for x in range(d)]
        #
        # nb = 2
        # t = 2
        # d = 10
        # w = 0.8
        # xp = [t * x + w * nb for x in range(d)]
        #
        # nb = 3
        # t = 3
        # d = 10
        # w = 0.8
        # xt = [t * x + w * nb for x in range(d)]
        # xlist = [xc, xp, xt]

    ## plot total value curve
    # ax = plt.subplot(121)
    # ax.plot(range(5), total_energ)
    # ax.plot(range(5), ela_energ)
    # ax.plot(range(5), pla_energ)
    # ax.legend(['total', 'elastic', 'plastic'])
    # ax.set_xticks(range(5))
    # ax.set_xticklabels(loadx)
    # ax.set(xlabel='loading (g)', ylabel='strain energy density(∆L/100gram)', title='Sequential strain energy density')

    ## plot incremental data
    # bx = plt.subplot(122)
    # bx = plt.subplot()
    # inc_t = [np.mean(t_inc_energ[loadx[1]]), np.mean(t_inc_energ[loadx[2]]), np.mean(t_inc_energ[loadx[3]]), np.mean(t_inc_energ[loadx[4]])]
    # inc_e = [np.mean(e_inc_energ[loadx[1]]), np.mean(e_inc_energ[loadx[2]]), np.mean(e_inc_energ[loadx[3]]), np.mean(e_inc_energ[loadx[4]])]
    # inc_p = [np.mean(p_inc_energ[loadx[1]]), np.mean(p_inc_energ[loadx[2]]), np.mean(p_inc_energ[loadx[3]]), np.mean(p_inc_energ[loadx[4]])]
    # er_t = [np.std(t_inc_energ[loadx[1]]), np.std(t_inc_energ[loadx[2]]), np.std(t_inc_energ[loadx[3]]), np.std(t_inc_energ[loadx[4]])]
    # er_e = [np.std(e_inc_energ[loadx[1]]), np.std(e_inc_energ[loadx[2]]), np.std(e_inc_energ[loadx[3]]), np.std(e_inc_energ[loadx[4]])]
    # er_p = [np.std(p_inc_energ[loadx[1]]), np.std(p_inc_energ[loadx[2]]), np.std(p_inc_energ[loadx[3]]), np.std(p_inc_energ[loadx[4]])]
    #
    # bx.bar(xt, inc_t, yerr = er_t)
    # bx.bar(xe, inc_e, yerr = er_e)
    # bx.bar(xp, inc_p, yerr = er_p)
    # bx.legend(['total', 'elastic', 'plastic'])
    # bx.set_xticks(xe)
    # bx.set_xticklabels(['8', '12', '16', '20'])
    # bx.set(xlabel='loading (g)', ylabel='Incremental strain energy density(∆L/100gram)', title='Incremental strain energy density')
    # plt.show()


    # output the results in a file
    # energ = pd.concat([t_energ, e_energ, p_energ], axis = 1)
    # inc_energ = pd.concat([t_inc_energ, e_inc_energ, p_inc_energ], axis = 1)
    # summary = pd.concat([energ, inc_energ], axis = 1)
    # summary = pd.concat([p_energ, hyst_energy, store_energ], axis = 1)
    # summary_inc = pd.concat([p_inc_energ, h_inc_energ, s_inc_energ], axis = 1)
    # print(summary)
    # summary.to_csv(output1)
    # summary_inc.to_csv(output2)

    C_p = []
    C_h = []
    C_s = []
    Ce_p = []
    Ce_h = []
    Ce_s = []
    for i in range(len(strainT)):
        C_p.append(np.mean(p_energ.iloc[:, i]))
        C_h.append(np.mean(hyst_energy.iloc[:, i]))
        C_s.append(np.mean(store_energ.iloc[:, i]))
        Ce_p.append(np.std(p_energ.iloc[:, i]))
        Ce_h.append(np.std(hyst_energy.iloc[:, i]))
        Ce_s.append(np.std(store_energ.iloc[:, i]))


    total_chart = pd.DataFrame(columns=strainT)
    total_chart.loc[len(total_chart)] = C_p
    total_chart.loc[len(total_chart)] = C_h
    total_chart.loc[len(total_chart)] = C_s

        # total_chart.to_csv(output)
    # Line-plot with error bar
    # ax = plt.subplot(121)
    # x = loadx
    # C_t = [np.mean(t_energ[loadx[0]]), np.mean(t_energ[loadx[1]]), np.mean(t_energ[loadx[2]]), np.mean(t_energ[loadx[3]]), np.mean(t_energ[loadx[4]])]
    # C_e = [np.mean(e_energ[loadx[0]]), np.mean(e_energ[loadx[1]]), np.mean(e_energ[loadx[2]]), np.mean(e_energ[loadx[3]]), np.mean(e_energ[loadx[4]])]
    # C_p = [np.mean(p_energ[loadx[0]]), np.mean(p_energ[loadx[1]]), np.mean(p_energ[loadx[2]]), np.mean(p_energ[loadx[3]]), np.mean(p_energ[loadx[4]])]
    # er_t = [np.std(t_energ[loadx[0]]), np.std(t_energ[loadx[1]]), np.std(t_energ[loadx[2]]), np.std(t_energ[loadx[3]]), np.std(t_energ[loadx[4]])]
    # er_e = [np.std(e_energ[loadx[0]]), np.std(e_energ[loadx[1]]), np.std(e_energ[loadx[2]]), np.std(e_energ[loadx[3]]), np.std(e_energ[loadx[4]])]
    # er_p = [np.std(p_energ[loadx[0]]), np.std(p_energ[loadx[1]]), np.std(p_energ[loadx[2]]), np.std(p_energ[loadx[3]]), np.std(p_energ[loadx[4]])]


    # plt.errorbar(x, C_t, yerr = er_t, color= 'blue')
    # plt.errorbar(x, C_e, yerr = er_e, color='green')
    # plt.errorbar(x, C_p, yerr = er_p, color='orange')
    #
    # ax.set_xticks(loadx)
    # ax.set_xticklabels(loadx)
    # ax.set(xlabel='loading (g)', ylabel='strain energy density(MPa)', title='Sequential strain energy density_' + name)
    # blue_line = mlines.Line2D([], [], color = 'blue', label = 'Input energy(total)' )
    # green_line = mlines.Line2D([], [], color = 'green', label = 'Stored energy(elasticity)' )
    # orange_line = mlines.Line2D([], [], color = 'orange', label = 'Dissipated energy(plasticity)' )
    # ax.legend(handles = [blue_line, green_line, orange_line], loc = 2)


    # plt.show()

    # bar-plot with error bar
    bx = plt.subplot(111)
    # x = [1, 2, 3, 4, 5, 6 , 7, 8, 9, 10]
    lw = 2
    ms = 6
    bx.set_xticks(strainT)
    bx.set_xticklabels(strainT)
    bx.errorbar(strainT, C_p, yerr = Ce_p, linestyle = ls[g], marker = '8', color='C9', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    bx.errorbar(strainT, C_h, yerr = Ce_h,  linestyle = ls[g], marker = '8', color='C3', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    bx.errorbar(strainT, C_s, yerr = Ce_s,  linestyle = ls[g], marker = '8', color='C2', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    blue_line = mlines.Line2D([], [], color = 'C2', label = 'Stored energy' )
    green_line = mlines.Line2D([], [], color = 'C3', label = 'Dissipated energy (elastic hysteresis)')
    orange_line = mlines.Line2D([], [], color = 'C9', label = 'Dissipated energy(plasticity)')
    bx.legend(handles = [blue_line, green_line, orange_line], loc = 2, fontsize=10, frameon = False)
    bx.set(xlabel='Strain (%)', ylabel='strain energy density(MPa)', title='Sequential strain energy density_' + name)


    # inc_t = [np.mean(t_inc_energ[loadx[1]]), np.mean(t_inc_energ[loadx[2]]), np.mean(t_inc_energ[loadx[3]]), np.mean(t_inc_energ[loadx[4]])]
    # inc_e = [np.mean(e_inc_energ[loadx[1]]), np.mean(e_inc_energ[loadx[2]]), np.mean(e_inc_energ[loadx[3]]), np.mean(e_inc_energ[loadx[4]])]
    # inc_p = [np.mean(p_inc_energ[loadx[1]]), np.mean(p_inc_energ[loadx[2]]), np.mean(p_inc_energ[loadx[3]]), np.mean(p_inc_energ[loadx[4]])]
    # er_t = [np.std(t_inc_energ[loadx[1]]), np.std(t_inc_energ[loadx[2]]), np.std(t_inc_energ[loadx[3]]), np.std(t_inc_energ[loadx[4]])]
    # er_e = [np.std(e_inc_energ[loadx[1]]), np.std(e_inc_energ[loadx[2]]), np.std(e_inc_energ[loadx[3]]), np.std(e_inc_energ[loadx[4]])]
    # er_p = [np.std(p_inc_energ[loadx[1]]), np.std(p_inc_energ[loadx[2]]), np.std(p_inc_energ[loadx[3]]), np.std(p_inc_energ[loadx[4]])]
    #
    # bx.bar(xt, inc_t, yerr = er_t, color = 'blue')
    # bx.bar(xe, inc_e, yerr = er_e, color = 'green')
    # bx.bar(xp, inc_p, yerr = er_p, color = 'orange')
    # bx.legend(['total', 'elastic', 'plastic'])
    # bx.set_xticks(xe)
    # bx.set_xticklabels(['8', '12', '16', '20'])
    # bx.set(xlabel='loading (g)', ylabel='Incremental strain energy density(MPa)', title='Incremental strain energy density_3onion')
    # plt.show()


    g += 1
# plt.axhline(y = 0, color = 'black')
plt.show()

