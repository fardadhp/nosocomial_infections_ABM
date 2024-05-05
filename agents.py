import pandas as pd
import numpy as np
import random
from collections import Counter
from datetime import datetime, timedelta
import os
import pickle

class Patient:
    def __init__(self, ID, hospital, unit, admissionDate, status, dischargeDate, contactPrecautions):
        self.ID = ID
        self.hospital = hospital
        self.unit = unit
        self.admissionDate = admissionDate
        self.status = status
        self.dischargeDate = dischargeDate
        self.LOS = (dischargeDate - admissionDate).days
        self.contactPrecautions = contactPrecautions
        self.treatmentStartDay = None
        self.room = None
        self.contactCount = 0
        self.active = True
        self.log = []

    def __del__(self):
        del self
    
    def removePatient(self, event):
        # save history
        if self.hospital.monteCarlo == False:
            pd.DataFrame(self.log, columns=['date','event']).to_csv(self.hospital.path+'/patients/patient_'+self.ID+'.csv', index=False)
        # remove all instances
        self.unit.patients.remove(self)
        self.unit.patients_list.remove(self.ID)
        self.unit = ''
        self.room = ''
        self.active = False
        if event == 'discharge':
            self.hospital.patients_list.remove(self.ID)
            self.hospital.patients.remove(self)
    
    def assignToRoom(self, adt_row):
        if type(adt_row['adt_room']) != type(str):
            self.room = Room(np.random.randint(100), self.unit, 0)
        else:    
            self.room = self.unit.rooms[adt_row['adt_room']]
            if type(self.room) == type(None):
                self.room = Room(adt_row['adt_room'], self.unit, 0)
        self.unit.rooms[adt_row['adt_room']] = self.room
        if self.hospital.output_verbose:
            self.log.append([adt_row['in_time'], f'assigned to room {self.room.ID}'])
   
    def infectionTreatment(self, date):
        if type(self.treatmentStartDay) == type(None):
            self.treatmentStartDay = date
            self.contactPrecautions = True
            if self.hospital.output_verbose:
                self.log.append([date, 'started antibiotics treatment'])
            # self.dischargeDate = max(self.dischargeDate, self.hospital.date+timedelta(days=self.unit.transmissionParams['infection_recovery']))
        # else:
        #     if date > self.treatmentStartDay + timedelta(days=self.unit.transmissionParams['infection_recovery']):
        #         self.treatmentStartDay = None
        #         if self.hospital.interventionType != 'universal_contact_precautions':
        #             self.contactPrecautions = False
        #         if self.hospital.output_verbose:
        #             self.log.append([date, 'ended antibiotics treatment'])
        #         if random.random() < self.unit.transmissionParams['proportion_of_infected_remain_colonized_after_treatment']:
        #             self.status = 'UC'
        #         else:
        #             self.status = 'X'

    def contactEnv(self, date):
        if self.hospital.transmissionPathways['patient_room']:
            # colonization (MOVED TO EVENT QUEUE)
            if self.status in ['S','X'] and self.room.contamination == 1:
                p = self.unit.transmissionParams['probability_environmental_colonization_hourly']
                if self.status == 'X':
                    p *= self.unit.transmissionParams['increase_factor_for_highly_susceptible']
                if random.random() < p:
                    self.status = 'UC'
                    if self.hospital.output_verbose:
                        self.log.append([date, f'colonized from environment at room {self.room.ID}'])
                    self.unit.log.append([date, self.ID, 'colonization', 'environment'])
            if self.status in ['UC','DC','I'] and self.room.contamination == 0:
                # shedding
                if random.random() < self.unit.transmissionParams['probability_room_contamination_by_colonized_patient']:
                    self.room.contamination = 1          

    def colonizationToInfection(self, date):
        if self.status in ['UC','DC']:
            time_of_infection = self.unit.transmissionParams['probability_colonization_to_infection'].sample(1)[0][0]
            if self.LOS > time_of_infection:
                self.status = 'I'
                if self.hospital.output_verbose:
                    self.log.append([date, 'infection symptom onset'])
                self.unit.log.append([date, self.ID, 'infection', 'colonization'])
                self.infectionTreatment(date)

    def backgroundTransmission(self, date):
        if random.random() < self.unit.transmissionParams['background_transmission_probability']:
            if random.random() < self.unit.transmissionParams['proportion_of_exposed_become_directly_infected']:
                self.status = 'I'
                if self.hospital.output_verbose:
                    self.log.append([date, 'infected from background'])
                self.unit.log.append([date, self.ID, 'infection', 'background'])
                self.infectionTreatment(date)
            else:
                self.status = 'UC'
                if self.hospital.output_verbose:
                    self.log.append([date, 'colonized from background'])
                self.unit.log.append([date, self.ID, 'colonization', 'background'])
                
    def surveillance(self, date):
        if (self.unit.surveillance_freq > 0) and ((self.hospital.date-self.hospital.T0).days % self.unit.surveillance_freq == 0):
            if self.status == 'UC':
                if random.random() < self.unit.prob_surveillance * self.unit.transmissionParams['test_sensitivity']:
                    self.status = 'DC'
                    if self.hospital.output_verbose:
                        self.log.append([date, 'colonization detected'])
                    self.unit.log.append([date, self.ID, 'detection', 'test'])
                    self.contactPrecautions = True
       
    def discharge(self, date, event, deceased=False):
        try:
            reason = f'discharged from {self.unit.name}'
        except:
            reason = 'discharged from unknown unit'
        if deceased:
            reason = 'deceased'
        if self.hospital.output_verbose:
            self.log.append([date, reason])
        if self.unit.terminal_room_disinfect_bin == 1:
            self.room.disinfect()
        self.removePatient(event)
        self.__del__()

        
class HCW:
    def __init__(self, ID, hospital, contamination, hygieneComplianceEnter, hygieneComplianceExit, PPECompliance):
        self.ID = ID
        self.hospital = hospital
        self.contamination = contamination
        self.hygieneComplianceEnter = hygieneComplianceEnter
        self.hygieneComplianceExit = hygieneComplianceExit
        self.PPECompliance = PPECompliance
        self.PPE = False
        self.PPEcontamination = 0
        self.log = []

    def washHands(self, hygieneCompliance, unit):
        # HCW hand hygiene
        if self.contamination > 0.99:
            if random.random() < hygieneCompliance:
                self.contamination = unit.transmissionParams['residual_contamination_post_hand_washing']
    
    def wearPPE(self):
        if random.random() < self.PPECompliance:
            self.PPE = True
            self.PPEcontamination = False
    
    def contactEnv(self, patient, date, env_to_hcw, hcw_to_env):
        if self.hospital.transmissionPathways['environmental']:
            # hcw contamination
            if patient.room.contamination == 1:
                if env_to_hcw:
                    if not self.PPE:
                        self.contamination = 1
                        if self.hospital.output_verbose:
                            self.log.append([date, f'contaminated by environment at unit {patient.unit.ID} room {patient.room.ID}'])
                    else:
                        self.PPEcontamination = 1
                        if self.hospital.output_verbose:
                            self.log.append([date, f'PPE contaminated by environment at unit {patient.unit.ID} room {patient.room.ID}'])
            # shedding
            else:
                if (not self.PPE and self.contamination == 1) or (self.PPE and self.PPEcontamination == 1):
                    if hcw_to_env:
                        patient.room.contamination = 1

    def interactWithPatient(self, patient, date, event_flags, p_hcw_to_patient):
        patient_to_hcw, hcw_to_patient, env_to_hcw, hcw_to_env = event_flags
        if patient.contactPrecautions == True:
            self.wearPPE()
        else:
            self.washHands(self.hygieneComplianceEnter, patient.unit)
        # how many times HCW contact env
        options = ['before & after', 'just before', 'just after', 'neither']
        probs = [0.25, 0.25, 0.25, 0.25]
        hcw_env_cont_times = options[np.where(np.random.multinomial(1, probs))[0][0]]
        # env <-> hcw at entry
        if 'before' in hcw_env_cont_times and any([env_to_hcw, hcw_to_env]):
            self.contactEnv(patient, date, env_to_hcw, hcw_to_env)
        # patient to HCW
        if patient.status in ['UC','DC','I']:
            if patient_to_hcw:
                if not self.PPE:
                    self.contamination = 1
                    if self.hospital.output_verbose:
                        self.log.append([date, f'contaminated by patient {patient.ID} at unit {patient.unit.ID} room {patient.room.ID}'])
                else:
                    self.PPEcontamination = 1
                    if self.hospital.output_verbose:
                        self.log.append([date, f'PPE contaminated by patient {patient.ID} at unit {patient.unit.ID} room {patient.room.ID}'])           
        # HCW to patient
        else:
            contamination = self.contamination
            if self.PPE:
                contamination = self.PPEcontamination
            probability_transmission = self.hospital.transmissionParams['probability_transmission_from_contaminated_hcw_to_susceptible_patient'] * contamination
            hcw_to_patient = (p_hcw_to_patient < probability_transmission)
            if (hcw_to_patient and (patient.status in ['S','X'])):
                patient.status = 'UC'
                if self.hospital.output_verbose:
                    patient.log.append([date, f'colonized by HCW {self.ID}'])
                    self.log.append([date, f'colonized patient {patient.ID}'])
                patient.unit.log.append([date, patient.ID, 'colonization', 'HCW'])
                        
        # env <-> hcw at exit
        if 'after' in hcw_env_cont_times and any([env_to_hcw, hcw_to_env]):
            self.contactEnv(patient, date, env_to_hcw, hcw_to_env)
        # remove PPE
        if self.PPE:
            self.PPE = False
            self.PPEcontamination = False
        else:
            self.washHands(self.hygieneComplianceExit, patient.unit)
            

    def contactStationBathroom(self, date):
        # shedding
        if self.contamination == 1:
            if self.hospital.transmissionPathways['nursing_station']:
                self.unit.station.contamination += self.unit.transmissionParams['proportion_of_exposed_become_directly_infected']
            if self.hospital.transmissionPathways['nurse_bathroom']:
                if random.random() < 1/6:
                    b = np.random.randint(len(self.unit.bathrooms))
                    self.unit.bathrooms[b].contamination += self.unit.transmissionParams['proportion_of_exposed_become_directly_infected']
        else:
            # contamination
            if self.hospital.transmissionPathways['nursing_station']:
                p = self.doseResponseFunction(self.unit.station.contamination)
                if random.random() < p:
                    self.contamination = 1
                    if self.hospital.output_verbose:
                        self.log.append([date, 'contaminated by environment at nursing station'])
            if self.hospital.transmissionPathways['nurse_bathroom']:
                if random.random() < 1/6:
                    b = np.random.randint(len(self.unit.bathrooms))
                    p = self.doseResponseFunction(self.unit.bathrooms[b].contamination)
                    if random.random() < p:
                        self.contamination = 1
                        if self.hospital.output_verbose:
                            self.log.append([date, 'contaminated by bathroom environment'])        


class Room:
    def __init__(self, ID, unit, contamination):
        self.ID = ID
        self.unit = unit
        self.contamination = contamination
        self.disinfect_efficacy = unit.transmissionParams['terminal_room_disinfection_efficacy']
        self.natural_clearance = unit.transmissionParams['pathogen_natural_clearance_from_dry_surfaces']
    
    def disinfect(self):
        if random.random() < self.disinfect_efficacy:
            self.contamination = 0
    
    def naturalClearance(self):
        if random.random() < self.natural_clearance:
            self.contamination = 0
    

class Unit:
    def __init__(self, hospital, hospital_ID, ID, name, capacity, bathrooms, \
        terminal_room_disinfect_bin, prob_adm_testing, prob_surveillance, \
            surveillance_freq, nurse_shift_hours):
        self.hospital = hospital
        self.hospital_ID = hospital_ID
        self.ID = ID
        self.name = name
        self.rooms = {}
        self.capacity = capacity
        self.bathrooms = bathrooms
        self.terminal_room_disinfect_bin = terminal_room_disinfect_bin
        self.prob_adm_testing = prob_adm_testing
        self.prob_surveillance = prob_surveillance
        self.surveillance_freq = surveillance_freq
        self.nurse_shift_hours = nurse_shift_hours
        self.transmissionParams = {}
        self.patients = {}
        self.patients_list = []
        self.station = None
        self.dailyAdmissionsCount = []
        self.dailyContacts = []
        self.stats = []
        self.log = []
        self.setup()

    def setup(self):
        # set transmission parameters
        self.transmissionParams = self.hospital.transmissionParams.copy()
        # create shared bathrooms
        b = self.bathrooms
        self.bathrooms = []
        for i in range(b):
            self.bathrooms.append(Room(i, self, 0))
        # create nursing station
        self.station = Room(0, self, 0)
        # create patient rooms
        rooms = self.hospital.event_queue.loc[self.hospital.event_queue['department_name_new']==self.name, 'adt_room'].dropna().unique()
        for r in rooms:
            self.rooms[r] = Room(r, self, np.random.randint(low=0, high=2))
        self.updateAdmissionDistribution()
        self.calculateComplianceBetaParameters()
       
    def makeNewPatient(self, adt_row):
        if adt_row['MaskedMRN'] in self.hospital.patients_list:
            admission_type = 'readmission'
        else:
            admission_type = 'admission'
        if self.hospital.randomImportation:
            inx = np.where(np.random.multinomial(1, self.transmissionParams['admission_distribution']))[0][0]
            status = 'SXCI'[inx]
            CP = False
            if status == "I":
                CP = True
            elif status == "C":
                if np.random.random() < self.prob_adm_testing * self.transmissionParams['test_sensitivity']:
                    status = 'DC'
                    CP = True
                else:
                    status = "UC"
        else:
            if any([adt_row[key] for key in ['MRSA_importation','VRE_importation']]):
                status = "DC"
                CP = True
            else:
                CP = False
                if random.random() < self.prob_adm_testing * (1 - self.transmissionParams['test_sensitivity']):
                    status = 'UC'
                elif random.random() < self.transmissionParams['highly_susceptible_ratio']:
                    status = 'X'
                else:
                    status = 'S'        
        if self.hospital.interventionType == "universal_contact_precautions":
            if adt_row['in_time'] >= self.hospital.interventionTime:
                CP = True
        if admission_type == 'admission':
            newPatient = Patient(adt_row['MaskedMRN'], self.hospital, self, adt_row['in_time'], status, adt_row['out_time'], CP)
        else:
            newPatient = self.hospital.patients[adt_row['MaskedMRN']]
            newPatient.unit = self
            newPatient.admissionDate = adt_row['in_time']
            newPatient.status = status
            newPatient.dischargeDate = adt_row['out_time']
            newPatient.LOS = (newPatient.dischargeDate - newPatient.admissionDate).days
            newPatient.contactPrecautions = CP 
            newPatient.active = True
        cases = {'S': [' as susceptible',''],
                 'X': [' as highly susceptible',''],
                 'UC': [' as undetected colonized','colonization'],
                 'DC': [' as detected colonized','colonization'],
                 'I': ['infected','infection']}
        if self.hospital.output_verbose:
            newPatient.log.append([adt_row['event_time'], f'{admission_type} as {cases[status][0]} to {self.name}'])
        if status in ['UC','DC','I']:                   
            newPatient.unit.log.append([adt_row['event_time'], adt_row['MaskedMRN'], cases[status][1], 'admission'])
        newPatient.assignToRoom(adt_row)
        self.patients_list.append(newPatient.ID)
        self.hospital.patients_list.append(newPatient.ID)
        self.patients[newPatient.ID] = newPatient
        self.hospital.patients[newPatient.ID] = newPatient
        return newPatient

    def calculateComplianceBetaParameters(self):
        # hygiene
        mu = self.transmissionParams["nurse_hygiene_compliance_enter_mean"]
        sigma = self.transmissionParams["HCW_hygiene_compliance_variability_from_mean"]
        mu_ = 1 / mu - 1
        sigma2 = sigma ** 2
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        self.transmissionParams['nurse_hygiene_compliance_enter_alpha'] = a
        self.transmissionParams['nurse_hygiene_compliance_enter_beta'] = a * mu_
        
        mu = self.transmissionParams["nurse_hygiene_compliance_exit_mean"]
        mu_ = 1 / mu - 1
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        self.transmissionParams['nurse_hygiene_compliance_exit_alpha'] = a
        self.transmissionParams['nurse_hygiene_compliance_exit_beta'] = a * mu_
        # PPE
        mu = self.transmissionParams["nurse_PPE_compliance_mean"]
        sigma = self.transmissionParams["HCW_PPE_compliance_variability_from_mean"]
        mu_ = 1 / mu - 1
        sigma2 = sigma ** 2
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        self.transmissionParams['nurse_PPE_compliance_alpha'] = a
        self.transmissionParams['nurse_PPE_compliance_beta'] = a * mu_
        
    def updateAdmissionDistribution(self):
        admission_C = self.transmissionParams['admission_prevalence']
        admission_I = 0
        admission_X = self.transmissionParams['highly_susceptible_ratio'] * (1 - admission_C - admission_I)
        admission_S = 1 - admission_C - admission_I - admission_X
        admission_S, admission_X, admission_C, admission_I = np.array([admission_S, admission_X, admission_C, admission_I]) / (admission_S + admission_X + admission_C + admission_I)
        self.transmissionParams['admission_distribution'] = [admission_S, admission_X, admission_C, admission_I]
    
    def replaceParameters(self, mc_params, rvs): 
        mc_params.loc[:len(rvs), 'value'] = rvs
        # replace sampled ICU parameters
        transm_mc_params = mc_params.loc[mc_params['owner']=='ICU',:].to_dict('records')
        for r in transm_mc_params:
            self.transmissionParams[r['model_param_name']] = r['value']    
        self.calculateComplianceBetaParameters()
        self.updateAdmissionDistribution()
    
    def findPatient(self, patient_ID):
        for patient in self.patients:
            if patient.ID == patient_ID:
                return patient
    
    def findRoom(self, room_ID):
        for room in self.rooms:
            if room.ID == room_ID:
                return room
    
    def discharge(self, adt_row):
        try:
            if adt_row['MaskedMRN'] in self.patients_list:
                patient = self.patients[adt_row['MaskedMRN']]
            else:
                patient = self.hospital.patients[adt_row['MaskedMRN']]
                patient.unit = self
                self.patients[adt_row['MaskedMRN']] = patient
                self.patients_list.append(patient.ID)
            room = self.rooms[adt_row['adt_room']]
            if patient.unit.name == self.name:
                patient.discharge(adt_row['event_time'], adt_row['event_type'])  
            room.disinfect() 
        except:
            pass

    def findPatient(self, patientID):
        for patient in self.patients:
            if patient.ID == patientID:
                return patient
                  
    def writeStats(self):
        lst = Counter([p.status for p in list(self.patients.values()) if p.active])
        lst = [lst[i] for i in ['S','X','UC','DC','I']]
        self.stats.append([*lst, self.dailyAdmissionsCount[-1], self.dailyContacts[-1]])

        
class Hospital:
    def __init__(self, ID, name, nUnits, simLength, burnIn, T0=None):
        self.ID = ID
        self.name = name
        self.nUnits = nUnits        
        self.simLength = simLength
        self.burnIn = burnIn
        self.T0 = T0
        self.units = {}
        self.HCWs = {}
        self.HCWs_list = []
        self.patients = {}
        self.patients_list = []
        self.cumNumPatients = 0
        self.transmissionParams = {}
        self.path = ''
        self.monteCarlo = False
        self.output_verbose = False
        self.fitToData = False
        self.randomImportation = True
        self.interventionType = 'None'
        self.interventionTime = None
        self.interventionEffect = 0
        self.eventcounter = 0
        self.setup()

    def calculateComplianceBetaParameters(self):
        # hygiene
        mu = self.transmissionParams["nurse_hygiene_compliance_enter_mean"]
        sigma = self.transmissionParams["HCW_hygiene_compliance_variability_from_mean"]
        mu_ = 1 / mu - 1
        sigma2 = sigma ** 2
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        if (a <= 0) | (mu_ <= 0):
            print("Beta:", mu, sigma)
        self.transmissionParams['nurse_hygiene_compliance_enter_alpha'] = a
        self.transmissionParams['nurse_hygiene_compliance_enter_beta'] = a * mu_
        
        mu = self.transmissionParams["nurse_hygiene_compliance_exit_mean"]
        mu_ = 1 / mu - 1
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        self.transmissionParams['nurse_hygiene_compliance_exit_alpha'] = a
        self.transmissionParams['nurse_hygiene_compliance_exit_beta'] = a * mu_
        # PPE
        mu = self.transmissionParams["nurse_PPE_compliance_mean"]
        sigma = self.transmissionParams["HCW_PPE_compliance_variability_from_mean"]
        mu_ = 1 / mu - 1
        sigma2 = sigma ** 2
        a = (mu_ - (1+mu_)**2 * sigma2) / ((1+mu_)**3 * sigma2)
        self.transmissionParams['nurse_PPE_compliance_alpha'] = a
        self.transmissionParams['nurse_PPE_compliance_beta'] = a * mu_
        
    def updateAdmissionDistribution(self):
        admission_C = self.transmissionParams['admission_prevalence']
        admission_I = 0
        admission_X = self.transmissionParams['highly_susceptible_ratio'] * (1 - admission_C - admission_I)
        admission_S = 1 - admission_C - admission_I - admission_X
        admission_S, admission_X, admission_C, admission_I = np.array([admission_S, admission_X, admission_C, admission_I]) / (admission_S + admission_X + admission_C + admission_I)
        self.transmissionParams['admission_distribution'] = [admission_S, admission_X, admission_C, admission_I]
        
    def setup(self):
        # transmission parameters
        params = pd.read_csv('./data/transmission_parameters.csv')
        paramsNames = params.loc[:,'parameters'].values
        paramsValues = params.loc[:,'value'].values
        for i,p in enumerate(paramsNames):
            self.transmissionParams[p] = paramsValues[i]
        with open("./data/kde_admission_to_infection_time_ALL.pkl", "rb") as f:
            self.transmissionParams['probability_colonization_to_infection'] = pickle.load(f)
        self.updateAdmissionDistribution()
        self.calculateComplianceBetaParameters()
        # transmission pathways
        self.transmissionPathways = dict(np.array(pd.read_csv('./data/transmission_pathways.csv', header=None)))
        # set T0
        if type(self.T0) == type(None):
            self.T0 = self.event_queue.loc[0, 'event_time']
        elif isinstance(self.T0, str):
            self.T0 = datetime.strptime(self.T0, '%Y-%m-%d')    
        self.date = self.T0 
        # import event queue data
        self.event_queue = pd.read_csv("./data/event_queue_all_2017_2018.csv")
        self.event_queue['MRSA_importation'] = self.event_queue['MRSA_importation'].fillna(0).astype(bool)
        self.event_queue['VRE_importation'] = self.event_queue['VRE_importation'].fillna(0).astype(bool) 
        for c in ['in_time','out_time','contact_time_start','contact_time_end','event_time']:
            self.event_queue[c] =  pd.to_datetime(self.event_queue[c], format='%Y-%m-%d  %H:%M:%S')       
        self.event_queue.loc[self.event_queue['event_time'].between(self.T0, self.T0+timedelta(days=self.simLength)), :]
        # self.addEnvColonizationEvents()
        # create units
        unitsData = pd.read_csv("./data/units_parameters.csv")        
        for i in range(unitsData.shape[0]):
            self.units[unitsData.loc[i, "name"]] = Unit(self, *unitsData.iloc[i,:].values)
        # set intervention(s)
        self.activateIntervention()
        # create output folders
        if len(self.path) == 0:
            self.path = './output/' + datetime.strftime(datetime.now(), '%Y_%m_%d_%H_%M')
        self.setOutputVerbose(self.output_verbose)
    
    # def addEnvColonizationEvents(self):
    #     eq = self.event_queue
    #     admissionEvents = eq.loc[eq['event_type']=="admission", :].to_dict('records')
    #     new_events = []
    #     for record in admissionEvents[:10]:
    #         los = int((record['out_time'] - record['in_time']).total_seconds() // 3600)
    #         if random.random() < 1 - (1 - self.transmissionParams['probability_environmental_colonization_hourly']) ** los:
    #             dtime = int(np.random.choice(np.arange(los)))
    #             time = record['in_time'] + timedelta(hours=dtime)
    #             record['event_type'] = 'env_colonization'
    #             record['event_time'] = time
    #             new_events.append(record)
    #     new_events = pd.DataFrame.from_dict(new_events)
    #     eq = pd.concat([eq, new_events])
    #     eq.sort_values(by=['event_time'], inplace=True)
    #     self.event_queue = eq
                
    def setOutputVerbose(self, output_verbose):
        self.output_verbose = output_verbose
        if output_verbose:
            paths = ['./output/', self.path, self.path+'/patients', self.path+'/units', self.path+'/HCWs']
            for p in paths:
                if not os.path.exists(p):
                    os.mkdir(p)
            
    
    def activateIntervention(self, intervention_name=''):
        intervention = pd.read_csv("./data/intervention_parameters.csv")
        if len(intervention_name) > 0:
            intervention = intervention.loc[intervention['type']==intervention_name, :]
        else:
            intervention = intervention.loc[intervention['activation']==1, :]
        if len(intervention) > 0:
            self.interventionType = intervention.type.values[0]
            self.interventionTime = intervention.time.values[0] + self.burnIn
            self.interventionEffect = intervention.impact.values[0]       
        if type(self.interventionTime) == type(""):
            self.interventionTime = datetime.strptime(self.interventionTime, '%Y-%m-%d') 
        elif type(self.interventionTime) == np.int64:
            self.interventionTime = self.T0 + timedelta(days=int(self.interventionTime))
        
    def decolonization(self, target_patients):
        for patient in target_patients:
            if random.random() < self.interventionEffect:
                patient.status = 'S'
                patient.contactPrecautions = False
    
    def reset(self):
        transmissionParams = self.transmissionParams
        units_transmissionParams = []
        for id, u in self.units.items():
            units_transmissionParams.append(u.transmissionParams)
        self.units = {}
        self.HCWs = {}
        self.HCWs_list = []
        self.patients = {}
        self.patients_list = []
        self.cumNumPatients = 0
        self.date = self.T0
        # recreate units
        unitsData = pd.read_csv("./data/units_parameters.csv")        
        for i in range(unitsData.shape[0]):
            self.units[unitsData.loc[i, "name"]] = Unit(self, *unitsData.iloc[i,:].values)
            self.units[unitsData.loc[i, "name"]].transmissionParams = units_transmissionParams[i]
            self.units[unitsData.loc[i, "name"]].updateAdmissionDistribution()
        self.updateAdmissionDistribution()
        self.transmissionParams = transmissionParams
            
    def findHCW(self, HCW_ID):
        for hcw in self.HCWs:
            if hcw.ID == HCW_ID:
                return hcw
    
    def findPatient(self, patient_ID):
        for patient in self.patients:
            if patient.ID == patient_ID:
                return patient
    
    def findUnit(self, unit_name):
        for unit in self.units:
            if unit.name == unit_name:
                return unit

    def startDay(self):        
        for id, unit in self.units.items():
            unit.dailyAdmissionsCount.append(0)
            # reset patients' daily counters + daily prob of developing an infection
            for id, patient in unit.patients.items():
                patient.contactCount = 0
                if patient.active:
                    patient.colonizationToInfection(self.date)
            # natural clearance
            for id, room in unit.rooms.items():
                room.naturalClearance()
            unit.station.naturalClearance()
            for bathroom in unit.bathrooms:
                bathroom.naturalClearance()
        # interventions
        if 'target_decolonization' in self.interventionType:
            # [print(p) for p in self.patients]
            colonized_patients = [self.patients[ID] for ID in self.patients_list if self.patients[ID].status == 'DC']
            if len(colonized_patients) > 0:
                self.decolonization(colonized_patients)

    def simulateDay(self):
        today_events = self.event_queue[((self.date<=self.event_queue['event_time'].values)&(self.event_queue['event_time'].values<self.date+timedelta(days=1)))]
        today_events.reset_index(drop=True, inplace=True)
        today_events = today_events.to_dict('records')
        # iterate over events
        for event in today_events:
            unit = self.units[event['department_name_new']]
            if event['event_type'] == 'admission':
                patient = unit.makeNewPatient(event)
                unit.dailyAdmissionsCount[-1] += 1 
            elif event['event_type'] in ['discharge','transfer']:
                unit.discharge(event)
            else:
                # contact
                event_probabilities = [unit.transmissionParams[c] for c in ['probability_HCW_contamination_from_colonized_patient',
                                                                            'probability_transmission_from_contaminated_hcw_to_susceptible_patient',
                                                                            'probability_hcw_contamination_from_contaminated_env',
                                                                            'probability_env_contamination_from_contaminated_hcw']]
                
                probs = np.random.random(len(event_probabilities))
                event_probabilities[0] *= unit.transmissionParams['shedding_increase_factor_for_infected']
                event_probabilities[1] *= unit.transmissionParams['increase_factor_for_highly_susceptible']
                event_flags = (probs < event_probabilities)
                # find patient
                if event['MaskedMRN'] in self.patients_list:
                    patient = self.patients[event['MaskedMRN']]
                else: # when patient's admission record not in the data
                    if type(event['department_name_new']) == type(''):
                        unit_name = event['department_name_new']
                    else:
                        unit_name = self.event_queue.loc[(self.event_queue['event_time']>=self.T0) &
                                                        (self.event_queue['event_type'].isin(['transfer','discharge'])) &
                                                        (self.event_queue['MaskedMRN']==event['MaskedMRN']),'department_name_new'].values
                        if len(unit_name) > 0:
                            unit_name = unit_name[0]
                    if len(unit_name) > 0:
                        my_unit = self.units[unit_name]
                    else:
                        my_unit = np.random.choice(self.units) 
                    patient = my_unit.makeNewPatient(event)
                    patient.unit = my_unit
                             
                if event['event_type'] == 'env_colonization':
                    patient.status = 'UC'
                    if self.output_verbose:
                        patient.log.append([event['event_time'], 'colonized from the room environment'])
                    patient.unit.log.append([event['event_time'], patient.ID, 'colonization', 'environment'])
                elif any(event_flags):
                    # find HCW
                    if event['HCW_ID'] in self.HCWs_list:
                        hcw = self.HCWs[event['HCW_ID']]
                    else:
                        alpha = unit.transmissionParams['nurse_hygiene_compliance_enter_alpha']
                        if alpha <= 0:
                            print(alpha)
                            print(unit.transmissionParams["nurse_hygiene_compliance_enter_mean"])
                            print(unit.transmissionParams["HCW_hygiene_compliance_variability_from_mean"])
                        nurse_hygiene_compliance_enter = np.random.beta(unit.transmissionParams['nurse_hygiene_compliance_enter_alpha'], unit.transmissionParams['nurse_hygiene_compliance_enter_beta'])
                        nurse_hygiene_compliance_exit = np.random.beta(unit.transmissionParams['nurse_hygiene_compliance_exit_alpha'], unit.transmissionParams['nurse_hygiene_compliance_exit_beta'])
                        nurse_PPE_compliance = np.random.beta(unit.transmissionParams['nurse_PPE_compliance_alpha'], unit.transmissionParams['nurse_PPE_compliance_beta'])
                        if "improve_low_compliance" in self.interventionType:
                            nurse_hygiene_compliance_enter = max(nurse_hygiene_compliance_enter, self.interventionEffect)
                            nurse_hygiene_compliance_exit = max(nurse_hygiene_compliance_exit, self.interventionEffect)
                            nurse_PPE_compliance = max(nurse_PPE_compliance, self.interventionEffect)
                        hcw = HCW(event['HCW_ID'], self, 0, nurse_hygiene_compliance_enter, nurse_hygiene_compliance_exit, nurse_PPE_compliance)
                        self.HCWs[event['HCW_ID']] = hcw
                        self.HCWs_list.append(event['HCW_ID'])
                    # HCW-patient interaction
                    if patient.active:
                        patient.contactCount += 1 
                        self.eventcounter += 1
                        hcw.interactWithPatient(patient, self.date, event_flags, probs[1])
                        # HCW transmission (prob[1]) has to be re-evaluated due to the residual contamination probability and compliance probability
                        # which vary by HCW (as we don't want to retireve the HCW agent to save time)
        # daily events 
        for id, patient in self.patients.items():
            if patient.active:
                patient.backgroundTransmission(self.date)
                patient.contactEnv(self.date)  
                patient.surveillance(self.date)     

    def endDay(self):
        for id, unit in self.units.items():  
            dailyContacts = 0
            for id, patient in unit.patients.items():
                dailyContacts += patient.contactCount
            unit.dailyContacts.append(dailyContacts)
            unit.writeStats()

    def simulate(self):
        for date in pd.date_range(start=self.T0, end=self.T0+timedelta(days=self.simLength)):
            self.date = date
            # print(f"New day: {date}")
            self.startDay()
            self.simulateDay()
            self.endDay()
        self.writeOutputs()
                                
    def writeOutputs(self):
        if (not self.monteCarlo) or self.output_verbose:
            # for patient in self.patients:
            #     pd.DataFrame(patient.log, columns=['date','event']).to_csv(self.path+'/patients/patient_ID_'+str(patient.ID)+'.csv', index=False)
            for u in self.units:
                pd.DataFrame(u.stats, columns=['S','X','UC','DC','I','admissions','contacts']).to_csv(self.path+'/units/unit_'+str(u.name)+'_stats.csv')
                pd.DataFrame(u.log, columns=['date','patient_ID','event','source']).to_csv(self.path+'/units/unit_'+str(u.name)+'_log.csv', index=False)
            # for hcw in self.HCWs:
            #     pd.DataFrame(hcw.log, columns=['date','event']).to_csv(self.path+'/HCWs/hcw_ID_'+str(hcw.ID).replace('/','')+'.csv', index=False)
            
