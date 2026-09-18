# physiology_builder.py
import numpy as np

class PhysiologyBuilder:
    def __init__(self):
        # Baseline parameters for a standard 70kg, healthy male
        self.base_weight = 70.0
        
        self.base_volumes = {
            'V_lungs': 3.3, 'V_adipose': 14.5, 'V_liver': 1.5,
            'V_muscle': 29.0, 'V_kidney': 0.3, 'V_skin': 2.6,
            'V_brain': 1.4, 'V_heart': 0.3, 'V_bone': 7.0,
            'V_blood': 5.0
        }
        
        self.base_flows = {
            'Q_lungs': 390.0, 'Q_adipose': 20.0, 'Q_liver': 96.0,
            'Q_muscle': 73.0, 'Q_kidney': 71.0, 'Q_skin': 22.0,
            'Q_brain': 46.0, 'Q_heart': 16.0, 'Q_bone': 19.0
        }

    def build_system(self, weight_kg=70.0, age=30, cirrhosis_status='None', renal_status='Healthy'):
        physio = {}
        
        # 1. Weight Scaling (Non-linear allometric scaling)
        weight_scalar = weight_kg / self.base_weight
        
        physio['V_adipose'] = self.base_volumes['V_adipose'] * (weight_scalar ** 1.2) # Obesity heavily expands fat
        physio['V_muscle'] = self.base_volumes['V_muscle'] * (weight_scalar ** 0.8)
        
        for tissue, vol in self.base_volumes.items():
            if tissue not in ['V_adipose', 'V_muscle']:
                physio[tissue] = vol * weight_scalar
                
        # Cardiac output and flows scale with weight^0.75
        co_scalar = weight_scalar ** 0.75
        for flow, q in self.base_flows.items():
            physio[flow] = q * co_scalar

        # 2. Disease State: Cirrhosis (Child-Pugh approximations)
        physio['Hepatic_CL_Scalar'] = 1.0
        if cirrhosis_status.lower() == 'mild':
            physio['Hepatic_CL_Scalar'] = 0.75
        elif cirrhosis_status.lower() == 'moderate':
            physio['Hepatic_CL_Scalar'] = 0.50
            physio['Q_liver'] *= 0.80 # Portal hypertension
        elif cirrhosis_status.lower() == 'severe':
            physio['Hepatic_CL_Scalar'] = 0.25
            physio['Q_liver'] *= 0.60
            physio['V_liver'] *= 0.70 # Liver shrinkage
            
        # 3. Disease State: Renal Impairment
        physio['Renal_CL_Scalar'] = 1.0
        if renal_status.lower() == 'ckd_stage_3':
            physio['Renal_CL_Scalar'] = 0.50
        elif renal_status.lower() == 'ckd_stage_4':
            physio['Renal_CL_Scalar'] = 0.25
            
        return physio