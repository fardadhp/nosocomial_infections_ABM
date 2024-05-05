import pandas as pd
import numpy as np
import random
from collections import Counter
from datetime import datetime, timedelta
import os
from tqdm import tqdm

class Patient:
    def __init__(self, ID, hospital, unit, admissionDate, status, dischargeDate, contactPrecautions):
        self.ID = ID
        self.hospital = hospital
        self.unit = unit
        self.admissionDate = admissionDate
        self.status = status
        self.dischargeDate = dischargeDate
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
            self.room = self.unit.findRoom(adt_row['adt_room'])
            if type(self.room) == type(None):
                self.room = Room(adt_row['adt_room'], self.unit, 0)
        self.unit.rooms.append(self.room)
        self.log.append([adt_row['in_time'], 'assigned to room '+str(self.room.ID)])
   
    def infectionTreatment(self, date):
        if type(self.treatmentStartDay) == type(None):
            self.treatmentStartDay = date
            self.contactPrecautions = True
            self.log.append([date, 'started antibiotics treatment'])
            self.dischargeDate = max(self.dischargeDate, self.hospital.date+timedelta(days=self.hospital.transmissionParams['infection_recovery']))
        else:
            if date > self.treatmentStartDay + timedelta(days=self.hospital.transmissionParams['infection_recovery']):
                self.treatmentStartDay = None
                if self.hospital.interventionType != 'universal_contact_precautions':
                    self.contactPrecautions = False
                self.log.append([date, 'ended antibiotics treatment'])
                if random.random() < self.hospital.transmissionParams['proportion_of_infected_remain_colonized_after_treatment']:
                    self.status = 'UC'
                else:
                    self.status = 'X'

    def contactEnv(self, date):
        if self.hospital.transmissionPathways['patient_room']:
            # colonization
            if self.status in ['S','X'] and self.room.contamination == 1:
                p = self.hospital.transmissionParams['probability_environmental_colonization']
                if self.status == 'X':
                    p *= self.hospital.transmissionParams['increase_factor_for_highly_susceptible']
                if random.random() < p:
                    self.status = 'UC'
                    self.log.append([date, 'colonized from environment at room '+str(self.room.ID)])
                    self.unit.log.append([date, self.ID, 'colonization', 'environment'])
            elif self.status in ['UC','DC'] and self.room.contamination == 0:
            # shedding
                if random.random() < self.hospital.transmissionParams['probability_room_contamination_by_colonized_patient']:
                    self.room.contamination = 1          

    def colonizationToInfection(self, date):
        if self.status in ['UC','DC']:
            if random.random() < self.hospital.transmissionParams['probability_colonization_to_infection']:
                self.status = 'I'
                self.log.append([date, 'infection symptom onset'])
                try:
                    self.unit.log.append([date, self.ID, 'infection', 'colonization'])
                except:
                    pass
                    # raise ValueError(self.ID, date)
                self.infectionTreatment(date)

    def backgroundTransmission(self, date):
        if random.random() < self.hospital.transmissionParams['background_transmission_probability']:
            if random.random() < self.hospital.transmissionParams['proportion_of_exposed_become_directly_infected']:
                self.status = 'I'
                self.log.append([date, 'infected from background'])
                self.unit.log.append([date, self.ID, 'infection', 'background'])
                self.infectionTreatment(date)
            else:
                self.status = 'UC'
                self.log.append([date, 'colonized from background'])
                self.unit.log.append([date, self.ID, 'colonization', 'background'])
                
    def surveillance(self, date):
        if (self.unit.surveillance_freq > 0) and ((self.hospital.date-self.hospital.T0).days % self.unit.surveillance_freq == 0):
            if self.status == 'UC':
                if random.random() < self.unit.prob_surveillance * self.hospital.transmissionParams['test_sensitivity']:
                    self.status = 'DC'
                    self.log.append([date, 'colonization detected'])
                    self.unit.log.append([date, self.ID, 'detection', 'test'])
                    self.contactPrecautions = True
       
    def discharge(self, date, event, deceased=False):
        try:
            reason = 'discharged from ' + self.unit.name
        except:
            pass
            # raise ValueError([self.ID, date])
        if deceased:
            reason = 'deceased'
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

    def washHands(self, hygieneCompliance):
        # HCW hand hygiene
        if self.contamination == 1:
            if random.random() < hygieneCompliance * self.hospital.transmissionParams['hands_hygiene_efficacy']:
                self.contamination = 0
    
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
                        self.log.append([date, 'contaminated by environment at unit '+str(patient.unit.ID)+' room '+str(patient.room.ID)])
                    else:
                        self.PPEcontamination = 1
                        self.log.append([date, 'PPE contaminated by environment at unit '+str(patient.unit.ID)+' room '+str(patient.room.ID)])
            # shedding
            else:
                if (not self.PPE and self.contamination == 1) or (self.PPE and self.PPEcontamination == 1):
                    if hcw_to_env:
                        patient.room.contamination = 1

    def interactWithPatient(self, patient, date, event_flags):
        patient_to_hcw, hcw_to_patient, env_to_hcw, hcw_to_env = event_flags
        if patient.contactPrecautions == 1:
            self.wearPPE()
        else:
            self.washHands(self.hygieneComplianceEnter)
        # env <-> hcw at entry
        if any([env_to_hcw, hcw_to_env]):
            self.contactEnv(patient, date, env_to_hcw, hcw_to_env)
        # patient to HCW
        if patient.status in ['UC','DC','I']:
            if patient_to_hcw:
                if not self.PPE:
                    self.contamination = 1
                    self.log.append([date, 'contaminated by patient '+str(patient.ID)+' at unit '+str(patient.unit.ID)+' room '+str(patient.room.ID)])
                else:
                    self.PPEcontamination = 1
                    self.log.append([date, 'PPE contaminated by patient '+str(patient.ID)+' at unit '+str(patient.unit.ID)+' room '+str(patient.room.ID)])           
        # HCW to patient
        else:
            contamination = self.contamination
            if self.PPE:
                contamination = self.PPEcontamination
            if patient.status in ['S','X'] and contamination == 1:
                if hcw_to_patient:
                    if random.random() < self.hospital.transmissionParams['proportion_of_exposed_become_directly_infected']:
                        patient.status = 'I'
                        patient.log.append([date, 'infected by HCW '+str(self.ID)])
                        patient.unit.log.append([date, patient.ID, 'infection', 'HCW'])
                        self.log.append([date, 'infected patient '+str(patient.ID)])
                        patient.infectionTreatment(date)
                    else:
                        patient.status = 'UC'
                        patient.log.append([date, 'colonized by HCW '+str(self.ID)])
                        patient.unit.log.append([date, patient.ID, 'colonization', 'HCW'])
                        self.log.append([date, 'colonized patient '+str(patient.ID)])
        # env <-> hcw at exit
        if any([env_to_hcw, hcw_to_env]):
            self.contactEnv(patient, date, env_to_hcw, hcw_to_env)
        # remove PPE
        if self.PPE:
            self.PPE = False
            self.PPEcontamination = False
        else:
            self.washHands(self.hygieneComplianceExit)
            

    def contactStationBathroom(self, date):
        # shedding
        if self.contamination == 1:
            if self.hospital.transmissionPathways['nursing_station']:
                self.unit.station.contamination += self.hospital.transmissionParams['proportion_of_exposed_become_directly_infected']
            if self.hospital.transmissionPathways['nurse_bathroom']:
                if random.random() < 1/6:
                    b = np.random.randint(len(self.unit.bathrooms))
                    self.unit.bathrooms[b].contamination += self.hospital.transmissionParams['proportion_of_exposed_become_directly_infected']
        else:
            # contamination
            if self.hospital.transmissionPathways['nursing_station']:
                p = self.doseResponseFunction(self.unit.station.contamination)
                if random.random() < p:
                    self.contamination = 1
                    self.log.append([date, 'contaminated by environment at nursing station'])
            if self.hospital.transmissionPathways['nurse_bathroom']:
                if random.random() < 1/6:
                    b = np.random.randint(len(self.unit.bathrooms))
                    p = self.doseResponseFunction(self.unit.bathrooms[b].contamination)
                    if random.random() < p:
                        self.contamination = 1
                        self.log.append([date, 'contaminated by bathroom environment'])        


class Room:
    def __init__(self, ID, unit, contamination):
        self.ID = ID
        self.unit = unit
        self.patients = []
        self.contamination = contamination
        self.disinfect_efficacy = unit.hospital.transmissionParams['disinfection_efficacy_for_dry_surfaces']
        self.natural_clearance = unit.hospital.transmissionParams['pathogen_natural_clearance_from_dry_surfaces']
    
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
        self.rooms = []
        self.capacity = capacity
        self.bathrooms = bathrooms
        self.terminal_room_disinfect_bin = terminal_room_disinfect_bin
        self.prob_adm_testing = prob_adm_testing
        self.prob_surveillance = prob_surveillance
        self.surveillance_freq = surveillance_freq
        self.nurse_shift_hours = nurse_shift_hours
        self.patients = []
        self.patients_list = []
        self.station = None
        self.dailyAdmissionsCount = []
        self.dailyContacts = []
        self.stats = []
        self.log = []
        self.setup()

    def setup(self):
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
            self.rooms.append(Room(r, self, np.random.randint(low=0, high=2)))
       
    def makeNewPatient(self, adt_row):
        if adt_row['MaskedMRN'] in self.hospital.patients_list:
            admission_type = 'readmission'
        else:
            admission_type = 'admission'
        if self.hospital.randomImportation:
            inx = np.where(np.random.multinomial(1, self.hospital.transmissionParams['admission_distribution']))[0][0]
            status = 'SXCI'[inx]
            CP = False
            if status == "I":
                CP = True
            elif status == "C":
                if np.random.random() < self.prob_adm_testing * self.hospital.transmissionParams['test_sensitivity']:
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
                if random.random() < self.prob_adm_testing * (1 - self.hospital.transmissionParams['test_sensitivity']):
                    status = 'UC'
                elif random.random() < self.hospital.transmissionParams['highly_susceptible_ratio']:
                    status = 'X'
                else:
                    status = 'S'        
        if self.hospital.interventionType == "universal_contact_precautions":
            if adt_row['in_time'] >= self.hospital.interventionTime:
                CP = True
        if admission_type == 'admission':
            newPatient = Patient(adt_row['MaskedMRN'], self.hospital, self, adt_row['in_time'], status, adt_row['out_time'], CP)
        else:
            newPatient = self.hospital.findPatient(adt_row['MaskedMRN'])   
            newPatient.unit = self
            newPatient.admissionDate = adt_row['in_time']
            newPatient.status = status
            newPatient.dischargeDate = adt_row['out_time']
            newPatient.contactPrecautions = CP 
            newPatient.active = True
        # if type(newPatient) == type(None):
        #     print(adt_row)
        #     raise TypeError("New patient was not created nor could be found.")
        cases = {'S': [' as susceptible',''],
                 'X': [' as highly susceptible',''],
                 'UC': [' as undetected colonized','colonization'],
                 'DC': [' as detected colonized','colonization'],
                 'I': ['infected','infection']}
        newPatient.log.append([adt_row['event_time'], admission_type+cases[status][0]+' to '+self.name])
        if status in ['UC','DC','I']:                   
            newPatient.unit.log.append([adt_row['event_time'], adt_row['MaskedMRN'], cases[status][1], 'admission'])
        newPatient.assignToRoom(adt_row)
        self.patients_list.append(newPatient.ID)
        self.hospital.patients_list.append(newPatient.ID)
        self.patients.append(newPatient)
        self.hospital.patients.append(newPatient)
        return newPatient
    
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
            patient = self.findPatient(adt_row['MaskedMRN'])
            if type(patient) == type(None):
                patient = self.hospital.findPatient(adt_row['MaskedMRN'])
                patient.unit = self
                self.patients.append(patient)
                self.patients_list.append(patient.ID)
            room = self.findRoom(adt_row['adt_room'])
            if patient.unit.name != self.name:
                pass
                # raise ValueError(["Patient's unit doesn't match discharge event record. Patient's log says': ", patient.unit.name, ", records says: ", self.name, ", for patient ", patient.ID])
            else:
                patient.discharge(adt_row['event_time'], adt_row['event_type'])  
        except:
            pass
        try:
            room.disinfect()         
        except:
            pass
            # raise ValueError(["Error in room disinfection proc.. records say room:", room, ", patient.room=", patient.room.ID, ", record: ", adt_row])

    def findPatient(self, patientID):
        for patient in self.patients:
            if patient.ID == patientID:
                return patient
      
    def writeStats(self):
        lst = Counter([p.status for p in self.patients if p.active])
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
        self.units = []
        self.HCWs = []
        self.HCWs_list = []
        self.patients = []
        self.patients_list = []
        self.cumNumPatients = 0
        self.transmissionParams = {}
        self.path = ''
        self.monteCarlo = False
        self.fitToData = False
        self.randomImportation = True
        self.interventionType = None
        self.interventionTime = None
        self.eventcounter = 0
        self.setup()

    def calculateUniformDistBoundaries(self, compliance_variability):
        self.transmissionParams['nurse_hygiene_compliance_enter_min'] = max(0, self.transmissionParams['nurse_hygiene_compliance_enter_mean'] - self.transmissionParams['nurse_hygiene_compliance_enter_half_interval'] * compliance_variability)
        self.transmissionParams['nurse_hygiene_compliance_enter_max'] = min(1, self.transmissionParams['nurse_hygiene_compliance_enter_mean'] + self.transmissionParams['nurse_hygiene_compliance_enter_half_interval'] * compliance_variability)
        self.transmissionParams['nurse_hygiene_compliance_exit_min'] = max(0, self.transmissionParams['nurse_hygiene_compliance_exit_mean'] - self.transmissionParams['nurse_hygiene_compliance_exit_half_interval'] * compliance_variability)
        self.transmissionParams['nurse_hygiene_compliance_exit_max'] = min(1, self.transmissionParams['nurse_hygiene_compliance_exit_mean'] + self.transmissionParams['nurse_hygiene_compliance_exit_half_interval'] * compliance_variability)
        self.transmissionParams['nurse_PPE_compliance_min'] = max(0, self.transmissionParams['nurse_PPE_compliance_mean'] - self.transmissionParams['nurse_PPE_compliance_half_interval'] * compliance_variability)
        self.transmissionParams['nurse_PPE_compliance_max'] = min(1, self.transmissionParams['nurse_PPE_compliance_mean'] + self.transmissionParams['nurse_PPE_compliance_half_interval'] * compliance_variability)       
        
    def updateAdmissionDistribution(self):
        admission_C = self.transmissionParams['admission_C']
        admission_I = self.transmissionParams['admission_I']
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
        self.updateAdmissionDistribution()
        self.calculateUniformDistBoundaries(1)
        # transmission pathways
        self.transmissionPathways = dict(np.array(pd.read_csv('./data/transmission_pathways.csv', header=None)))
        # import event queue data
        self.event_queue = pd.read_csv("./data/event_queue_all_2017_2018.csv")
        self.event_queue['MRSA_importation'] = self.event_queue['MRSA_importation'].fillna(0).astype(bool)
        self.event_queue['VRE_importation'] = self.event_queue['VRE_importation'].fillna(0).astype(bool) 
        for c in ['in_time','out_time','contact_time_start','contact_time_end','event_time']:
            self.event_queue[c] =  pd.to_datetime(self.event_queue[c], format='%Y-%m-%d  %H:%M:%S')       
        # create units
        unitsData = pd.read_csv("./data/units_parameters.csv")        
        for i in range(unitsData.shape[0]):
            self.units.append(Unit(self, *unitsData.iloc[i,:].values)) 
        # set T0
        if type(self.T0) == type(None):
            self.T0 = self.event_queue.loc[0, 'event_time']
        elif type(self.T0) == type(''):
            self.T0 = datetime.strptime(self.T0, '%Y-%m-%d')    
        self.date = self.T0 
        # set intervention(s)
        self.activateIntervention()
        # create output folders
        if len(self.path) == 0:
            self.path = './output/' + datetime.strftime(datetime.now(), '%Y_%m_%d_%H_%M')
        try:
            os.mkdir('./output/')
            os.mkdir(self.path)
            os.mkdir(self.path+'/patients')
            os.mkdir(self.path+'/units')
            os.mkdir(self.path+'/HCWs')
        except:
            pass
    
    def activateIntervention(self):
        intervention = pd.read_csv("./data/intervention_parameters.csv")
        intervention = intervention.loc[intervention['activation']==1, :]
        if len(intervention) > 0:
            self.interventionType = intervention.type.values[0]
            self.interventionTime = intervention.time.values[0] + self.burnIn           
        if type(self.interventionTime) == type(''):
            self.interventionTime = datetime.strptime(self.interventionTime, '%Y-%m-%d') 
        elif type(self.interventionTime) == np.int64:
            self.interventionTime = self.T0 + timedelta(days=int(self.interventionTime))
        
    def reset(self):
        self.units = []
        self.HCWs = []
        self.HCWs_list = []
        self.patients = []
        self.patients_list = []
        self.cumNumPatients = 0
        self.date = self.T0
        # recreate units
        unitsData = pd.read_csv("./data/units_parameters.csv")        
        for i in range(unitsData.shape[0]):
            self.units.append(Unit(self, *unitsData.iloc[i,:].values)) 
        self.updateAdmissionDistribution()
            
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
        for unit in self.units:
            unit.dailyAdmissionsCount.append(0)
            # reset patients' daily counters + daily prob of developing an infection
            for patient in unit.patients:
                patient.contactCount = 0
                patient.colonizationToInfection(self.date)
            # natural clearance
            for room in unit.rooms:
                room.naturalClearance()
            unit.station.naturalClearance()
            for bathroom in unit.bathrooms:
                bathroom.naturalClearance()

    def simulateDay(self):
        today_events = self.event_queue[((self.date<=self.event_queue['event_time'].values)&(self.event_queue['event_time'].values<self.date+timedelta(days=1)))]
        today_events.reset_index(drop=True, inplace=True)
        today_events = today_events.to_dict('records')
        # iterate over events
        for event in today_events:
            unit = self.findUnit(event['department_name_new'])
            if event['event_type'] == 'admission':
                patient = unit.makeNewPatient(event)
                unit.dailyAdmissionsCount[-1] += 1 
            elif event['event_type'] in ['discharge','transfer']:
                unit.discharge(event)
            else:
                # contact
                event_probabilities = [self.transmissionParams[c] for c in ['probability_transmission_from_infected_patient_to_uncontaminated_hcw',
                                                                            'probability_transmission_from_contaminated_hcw_to_susceptible_patient',
                                                                            'probability_hcw_contamination_from_contaminated_env',
                                                                            'probability_env_contamination_from_contaminated_hcw']]
                
                probs = np.random.random(len(event_probabilities))
                event_probabilities[0] *= self.transmissionParams['shedding_increase_factor_for_infected']
                event_probabilities[1] *= self.transmissionParams['increase_factor_for_highly_susceptible']
                event_flags = (probs < event_probabilities)
                if any(event_flags):
                    self.eventcounter += 1
                    # find patient
                    patient = self.findPatient(event['MaskedMRN'])
                    if type(patient) == type(None): # when patient's admission record not in the data
                        if type(event['department_name_new']) == type(str):
                            unit_name = event['department_name_new']
                        else:
                            unit_name = self.event_queue.loc[(self.event_queue['event_time']>=self.T0) &
                                                            (self.event_queue['event_type'].isin(['transfer','discharge'])) &
                                                            (self.event_queue['MaskedMRN']==event['MaskedMRN']),'department_name_new'].values
                            if len(unit_name) > 0:
                                unit_name = unit_name[0]
                        if len(unit_name) > 0:
                            for unit in self.units:
                                if unit.name == unit_name:
                                    my_unit = unit
                                    break 
                        else:
                            my_unit = np.random.choice(self.units) 
                        try:
                            patient = my_unit.makeNewPatient(event)
                        except:
                            pass
                            # raise ValueError(event)
                    # determine if any transmission events may happen 
                    # (added to improve performance but has marginal effect cause of other influencing parameters such as the status of patients and HCWs)
                    patient.contactCount += 1              
                
                
                
                    # find HCW
                    hcw = self.findHCW(event['HCW_ID'])
                    if type(hcw) == type(None):
                        nurse_hygiene_compliance_enter = np.random.uniform(self.transmissionParams['nurse_hygiene_compliance_enter_min'], self.transmissionParams['nurse_hygiene_compliance_enter_max'])
                        nurse_hygiene_compliance_exit = np.random.uniform(self.transmissionParams['nurse_hygiene_compliance_exit_min'], self.transmissionParams['nurse_hygiene_compliance_exit_max'])
                        nurse_PPE_compliance = np.random.uniform(self.transmissionParams['nurse_PPE_compliance_min'], self.transmissionParams['nurse_PPE_compliance_max'])
                        hcw = HCW(event['HCW_ID'], self, 0, nurse_hygiene_compliance_enter, nurse_hygiene_compliance_exit, nurse_PPE_compliance)
                        self.HCWs.append(hcw)
                        self.HCWs_list.append(event['HCW_ID'])
                    # HCW-patient interaction
                    if patient.active:
                        hcw.interactWithPatient(patient, self.date, event_flags)
                    else:
                        pass
                        # raise ValueError(["Contact with inactive patient.", event])
        # daily events 
        for patient in self.patients:
            if patient.active:
                patient.backgroundTransmission(self.date)
                patient.contactEnv(self.date)  
                patient.surveillance(self.date)     
                ## end infection treatment
                if patient.status == 'I':
                    patient.infectionTreatment(self.date)

    def endDay(self):
        for unit in self.units:  
            dailyContacts = 0
            for patient in unit.patients:
                dailyContacts += patient.contactCount
            unit.dailyContacts.append(dailyContacts)
            unit.writeStats()

    def simulate(self):
        if self.monteCarlo == False:
            try:
                os.mkdir(self.path)
                os.mkdir(self.path+'/patients')
                os.mkdir(self.path+'/units')
                os.mkdir(self.path+'/HCWs')
            except:
                pass
        for date in pd.date_range(start=self.T0, end=self.T0+timedelta(days=self.simLength)):
            self.date = date
            self.startDay()
            self.simulateDay()
            self.endDay()
        self.writeOutputs()
    
    def readFittingData(self, pathogens):
        for unit in self.units:
            data = []
            for pathogen in pathogens:
                data.append(pd.read_csv("./calibration/data/"+pathogen+"_acquisition_"+unit.name+".csv"))
            if len(data) > 1:
                data = pd.concat(data)
                data = data.fillna(0).reset_index(drop=True)
                data = data.groupby(by=['Unnamed: 0']).sum()
                data['total_acquisitions'] = data.sum(1)
                unit.data = data[['total_acquisitions']]
                                
    def writeOutputs(self):
        if self.monteCarlo == False:
            for patient in self.patients:
                pd.DataFrame(patient.log, columns=['date','event']).to_csv(self.path+'/patients/patient_ID_'+str(patient.ID)+'.csv', index=False)
            for u in self.units:
                pd.DataFrame(u.stats, columns=['S','X','UC','DC','I','admissions','contacts']).to_csv(self.path+'/units/unit_'+str(u.name)+'_stats.csv')
                pd.DataFrame(u.log, columns=['date','patient_ID','event','source']).to_csv(self.path+'/units/unit_'+str(u.name)+'_log.csv', index=False)
            for hcw in self.HCWs:
                pd.DataFrame(hcw.log, columns=['date','event']).to_csv(self.path+'/HCWs/hcw_ID_'+str(hcw.ID).replace('/','')+'.csv', index=False)
            
