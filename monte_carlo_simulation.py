import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
from tqdm import tqdm
from agents import Hospital
from datetime import datetime, timedelta
from scipy.stats import qmc
from glob import glob
import pickle
import argparse

def runMonteCarlo(inputList):
    samples_split, mc_params, path, iterations, days, burnIn, T0, pn = inputList
    hospitalsData = pd.read_csv("./data/hospitals_data.csv")
    nUnits = hospitalsData['number_units'].values[0]
    S, X, UC, DC, I, daily_admissions, daily_contacts, env, hcw, importation, admC, admI, \
    incC, incI, UA, qcol_per_1000_patient_days_pre_intervention, \
    qcol_per_1000_patients_pre_intervention, qcol_per_1000_contacts_pre_intervention, \
    qcol_per_1000_patient_days_post_intervention, qcol_per_1000_patients_post_intervention, \
    qcol_per_1000_contacts_post_intervention, qcol_acq_count_pre_intervention, qcol_acq_count_post_intervention, \
    qinf_per_1000_patient_days_pre_intervention, qinf_per_1000_patients_pre_intervention, qinf_per_1000_contacts_pre_intervention, \
    qinf_per_1000_patient_days_post_intervention, qinf_per_1000_patients_post_intervention, \
    qinf_per_1000_contacts_post_intervention, qinf_acq_count_pre_intervention, qinf_acq_count_post_intervention = [[[] for i in range(nUnits)] for j in range(31)]
    for rvs in tqdm(samples_split):
        mc_params['value'] = rvs
        hospital = Hospital(*hospitalsData.iloc[0,:].values, days, burnIn, T0)
        hospital.monteCarlo = True
        hospital.path = path
        hospital.monteCarloResults = {}
        hospital.randomImportation = True
        # replace sampled ICU parameters
        transm_mc_params = mc_params.loc[mc_params['owner']=='ICU',:].to_dict('records')
        for r in transm_mc_params:
            hospital.transmissionParams[r['model_param_name']] = r['value']    
        # replace other parameters (manually)
        if 'HCW compliance variability from mean' in mc_params['parameter']:
            compliance_variability = mc_params.loc[mc_params['parameter']=='HCW compliance variability from mean','value'].values[0]
            hospital.calculateUniformDistBoundaries(compliance_variability)
        if 'probability of contamination' in mc_params['parameter'].values:
            contamination_params = [
                'probability_room_contamination_by_colonized_patient',
                'probability_transmission_from_infected_patient_to_uncontaminated_hcw',
                'probability_hcw_contamination_from_contaminated_env',
                'probability_env_contamination_from_contaminated_hcw'
            ]
            ind1 = np.where(mc_params['model_param_name'].isin(contamination_params))
            ind2 = np.where(mc_params['model_param_name']=='probability_of_contamination')[0][0]
            rvs[ind1] = rvs[ind2]
            rvs = np.delete(rvs, ind2) 
            for p in contamination_params:
                hospital.transmissionParams[p] = mc_params.loc[mc_params['parameter']=='probability of contamination','value'].values[0]
        for i in range(iterations):
            if i > 0:
                hospital.reset()
            hospital.simulate()
            dates = pd.date_range(hospital.T0.date(), hospital.date.date()).date
            for j, u in enumerate(hospital.units):
                S[j].append(np.array(u.stats)[:,0])
                X[j].append(np.array(u.stats)[:,1])
                UC[j].append(np.array(u.stats)[:,2])
                DC[j].append(np.array(u.stats)[:,3])
                I[j].append(np.array(u.stats)[:,4])
                daily_admissions[j].append(np.array(u.stats)[:,5])
                daily_contacts[j].append(np.array(u.stats)[:,6])                
                # calculate contribution
                data = pd.DataFrame(u.log, columns=['date','patient_ID','event','source'])
                # data['date'] = pd.to_datetime(data['date']).apply(lambda x: x.date())
                data['date'] = pd.to_datetime(data['date']).dt.date
                data = data.loc[data['date']>=hospital.T0.date()+timedelta(hospital.burnIn),:].reset_index(drop=True)
                data = data.loc[data['event'].isin(['infection','colonization']),:]
                data.reset_index(drop=True, inplace=True)
                contribution = []
                for pathway in ['environment','HCW','admission']:
                    count = data.loc[data['source']==pathway,:].shape[0]
                    contribution.append(count)
                nEvents = max(data.shape[0], 1)
                contribution = np.array(contribution) / nEvents * 100
                env[j].append(contribution[0])
                hcw[j].append(contribution[1])
                importation[j].append(contribution[2])
                # calculate incidence nad importation
                data = pd.DataFrame(u.log, columns=['date','patient_ID','event','source'])
                # data['date'] = pd.to_datetime(data['date']).apply(lambda x: x.date())
                data['date'] = pd.to_datetime(data['date']).dt.date
                data = data.loc[data['date']>=(hospital.T0+timedelta(hospital.burnIn)).date(),:].reset_index(drop=True)              
                cols = ['colonization_importation', 'infection_importation', \
                'colonization_incidence', 'infection_incidence']
                incidence = pd.DataFrame(np.zeros((hospital.simLength-hospital.burnIn, len(cols))), index=dates[hospital.burnIn+1:], columns=cols)
                for day in incidence.index:
                    subset = data[((day<=data['date'].values)&(data['date'].values<day+timedelta(days=1)))]
                    if len(subset) > 0:
                        incidence.loc[day,'colonization_importation'] = subset.loc[(subset['event']=='colonization')&(subset['source']=='admission')].shape[0]
                        incidence.loc[day,'infection_importation'] = subset.loc[(subset['event']=='infection')&(subset['source']=='admission')].shape[0]
                        incidence.loc[day,'colonization_incidence'] = subset.loc[(subset['event']=='colonization')&(subset['source'].isin(['HCW','environment']))].shape[0]
                        incidence.loc[day,'infection_incidence'] = subset.loc[(subset['event']=='infection')&(subset['source'].isin(['HCW','environment','colonization']))].shape[0]
                admC[j].append(incidence.iloc[:,0].sum())
                admI[j].append(incidence.iloc[:,1].sum())
                incC[j].append(incidence.iloc[:,2].sum())
                incI[j].append(incidence.iloc[:,3].sum())
                # incidence rate
                stats = pd.DataFrame(u.stats, columns=['S','X','UC','DC','I','admissions','contacts'])
                census = stats[['S','X','UC','DC','I']].iloc[burnIn:,:].sum(1).values
                hospitalization = stats[['admissions']].iloc[burnIn:,0].values
                contacts = stats[['contacts']].iloc[burnIn:,0].values
                r = [[] for i in range(16)]
                if type(hospital.interventionTime) == type(None):
                    intv_time_index = 1
                else:
                    intv_time_index = (hospital.interventionTime - hospital.T0).days - hospital.burnIn
                with np.errstate(all='ignore'):
                    r[3].append(incidence['colonization_incidence'].values[intv_time_index:].sum() / census[intv_time_index:].sum() * 1000)
                    r[4].append(incidence['colonization_incidence'].values[intv_time_index:].sum() / hospitalization[intv_time_index:].sum() * 1000)
                    r[5].append(incidence['colonization_incidence'].values[intv_time_index:].sum() / contacts[intv_time_index:].sum() * 1000)                
                    r[7].append(incidence['colonization_incidence'].values[intv_time_index:].sum())
                    r[11].append(incidence['infection_incidence'].values[intv_time_index:].sum() / census[intv_time_index:].sum() * 1000)
                    r[12].append(incidence['infection_incidence'].values[intv_time_index:].sum() / hospitalization[intv_time_index:].sum() * 1000)
                    r[13].append(incidence['infection_incidence'].values[intv_time_index:].sum() / contacts[intv_time_index:].sum() * 1000)               
                    r[15].append(incidence['infection_incidence'].values[intv_time_index:].sum())
                if type(hospital.interventionTime) == type(None):
                    r[0] = r[3]
                    r[1] = r[4]
                    r[2] = r[5] 
                    r[6] = r[7]
                    r[8] = r[11]
                    r[9] = r[12]
                    r[10] = r[13] 
                    r[14] = r[15]
                else:
                    with np.errstate(all='ignore'):
                        r[0].append(incidence['colonization_incidence'].values[:intv_time_index].sum() / census[:intv_time_index].sum() * 1000)
                        r[1].append(incidence['colonization_incidence'].values[:intv_time_index].sum() / hospitalization[:intv_time_index].sum() * 1000)
                        r[2].append(incidence['colonization_incidence'].values[:intv_time_index].sum() / contacts[:intv_time_index].sum() * 1000)
                        r[6].append(incidence['colonization_incidence'].values[:intv_time_index].sum())
                        r[8].append(incidence['infection_incidence'].values[:intv_time_index].sum() / census[:intv_time_index].sum() * 1000)
                        r[9].append(incidence['infection_incidence'].values[:intv_time_index].sum() / hospitalization[:intv_time_index].sum() * 1000)
                        r[10].append(incidence['infection_incidence'].values[:intv_time_index].sum() / contacts[:intv_time_index].sum() * 1000)
                        r[14].append(incidence['infection_incidence'].values[:intv_time_index].sum())               
                r = np.array(r)
                r[np.isnan(r)] = 0
                UA[j].append(np.hstack([rvs, *r]))
    
    # write to file
    for i in range(nUnits):
        pd.DataFrame(np.array(S[i]).T).to_csv(path+'/S_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(X[i]).T).to_csv(path+'/X_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(UC[i]).T).to_csv(path+'/UC_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(DC[i]).T).to_csv(path+'/DC_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(I[i]).T).to_csv(path+'/I_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(daily_admissions[i]).T).to_csv(path+'/daily_admissions_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(daily_contacts[i]).T).to_csv(path+'/daily_contacts_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(env[i]).T).to_csv(path+'/env_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(hcw[i]).T).to_csv(path+'/hcw_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(importation[i]).T).to_csv(path+'/import_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(admC[i]).T).to_csv(path+'/admC_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(admI[i]).T).to_csv(path+'/admI_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(incC[i]).T).to_csv(path+'/incC_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(incI[i]).T).to_csv(path+'/incI_'+str(i)+'_'+str(pn)+'.csv', index=False)
        pd.DataFrame(np.array(UA[i]).transpose()).to_csv(path+'/UA_'+str(i)+'_'+str(pn)+'.csv', index=False)


def sampledWithIteration(nsamples, iterations, T0, burnIn, simLength, calibration=False):
    numproc = mp.cpu_count()
    days = simLength + burnIn
    if calibration:
        mc_range = pd.read_csv("./calibration/mc_parameters_ranges.csv")
        path = './calibration/' + datetime.strftime(datetime.now(), '%Y_%m_%d_%H_%M')
    else:
        mc_range = pd.read_csv("./data/mc_parameters_ranges.csv")
        path = './monte carlo/' + datetime.strftime(datetime.now(), '%Y_%m_%d_%H_%M')
    mc_params = mc_range.loc[:,['parameter','owner','model_param_name']]
    try:
        os.mkdir(path)
    except:
        pass
    lower_bounds = mc_range.loc[:,'min'].values
    upper_bounds = mc_range.loc[:,'max'].values
    sampler = qmc.LatinHypercube(d=len(mc_range))
    samples = sampler.random(n=nsamples)
    samples = qmc.scale(samples, lower_bounds, upper_bounds)
    if nsamples > 1:
        samples_split = np.array_split(samples, numproc)
    else:
        samples_split = [samples for i in range(numproc)]
        iterations  = int(iterations / numproc) + 1
    inputList = []
    for i in range(len(samples_split)):
        inputList.append([samples_split[i], mc_params, path, iterations, days, burnIn, T0, i])
    pools = mp.Pool(processes=numproc)
    pools.map(runMonteCarlo, inputList)
    # aggregate
    mc_params = mc_params.drop(mc_params[mc_params['parameter']=='probability of contamination'].index)
    mc_params.reset_index(drop=True, inplace=True)
    cols = ['S','X','UC','DC','I', 'daily_admissions', 'daily_contacts','env','hcw','import','admC', \
            'admI','incC','incI']#, 'qcol_rpd_pre','qcol_rp_pre','qcol_rc_pre', \
            # 'qcol_rpd_post','qcol_rp_post','qcol_rc_post','qcol_count_pre','qcol_count_post', \
            # 'qinf_rpd_pre','qinf_rp_pre','qinf_rc_pre', \
            # 'qinf_rpd_post','qinf_rp_post','qinf_rc_post','qinf_count_pre','qinf_count_post']
    monteCarloResults = {}
    nUnits = int(len(glob(path+'/admC*')) / numproc)
    for c in cols:
        l = [[] for u in range(nUnits)]
        for u in range(nUnits):
            for pn in range(numproc):
                filename = path+'/'+c+'_'+str(u)+'_'+str(pn)+'.csv'
                l[u].extend(pd.read_csv(filename).T.values)
                os.remove(filename)
            pd.DataFrame(l[u]).T.to_csv(path+'/'+c+'_'+str(u)+'.csv', index=False)
        monteCarloResults[c] = np.array(l)
    # uncertainty analysis
    UA = [[] for u in range(nUnits)]
    for u in range(nUnits):
        UA[u] = mc_params.copy()
        for c in ['qcol_per_1000_patient_days_pre_intervention',
                  'qcol_per_1000_patients_pre_intervention',
                  'qcol_per_1000_contacts_pre_intervention',
                  'qcol_per_1000_patient_days_post_intervention',
                  'qcol_per_1000_patients_post_intervention',
                  'qcol_per_1000_contacts_post_intervention',
                  'qcol_acq_count_pre_intervention',
                  'qcol_acq_count_post_intervention',
                  'qinf_per_1000_patient_days_pre_intervention',
                  'qinf_per_1000_patients_pre_intervention',
                  'qinf_per_1000_contacts_pre_intervention',
                  'qinf_per_1000_patient_days_post_intervention',
                  'qinf_per_1000_patients_post_intervention',
                  'qinf_per_1000_contacts_post_intervention',
                  'qinf_acq_count_pre_intervention',
                  'qinf_acq_count_post_intervention']:
            UA[u].loc[len(UA[u])] = [c, 'ICU', c]
        UA[u] = UA[u].T.values.tolist()
        for pn in range(numproc):
            filename = path+'/UA_'+str(u)+'_'+str(pn)+'.csv'
            UA[u].extend(pd.read_csv(filename).T.values)
            os.remove(filename)
        pd.DataFrame(UA[u]).T.to_csv(path+'/UA_'+str(u)+'.csv', index=False)
    monteCarloResults['UA'] = UA
    monteCarloResults['T0'] = T0
    monteCarloResults['burnIn'] = burnIn
    monteCarloResults['days'] = days
    monteCarloResults['nsamples'] = nsamples
    monteCarloResults['iterations'] = iterations
    
    with open(path+'/monteCarloResults.pickle', 'wb') as handle:
        pickle.dump(monteCarloResults, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    try:
        os.mkdir('monte carlo')
    except:
        pass
    parser = argparse.ArgumentParser(description="Help",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--nsamples", required=True, type=int, help="integer, number of set of samples (draws from distributions)", default=100)
    parser.add_argument("-i", "--iterations", type=int, help="integer, number of simulations per set of samples", default=30)
    parser.add_argument("-t", "--T0", type=str, help="string, starting date of simulation, format: 'YYYY-MM-DD'", default=None)
    parser.add_argument("-b", "--burnIn", type=int, help="integer, burn-in period", default=30)
    parser.add_argument("-l", "--simLength", type=int, help="integer, simulation length", default=360)
    parser.add_argument("-c", "--ifCalibration", type=bool, help="True/False", default=True)
    args = parser.parse_args()
    config = vars(args)
    nsamples, iterations, T0, burnIn, simLength, ifCalibration = config.values()

    sampledWithIteration(nsamples, iterations, T0, burnIn, simLength, calibration=ifCalibration)
