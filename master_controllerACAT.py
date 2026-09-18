# master_controllerACAT.py
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import argparse
import sys
import os

from distribution import calculate_tissue_kp
from physiology_builder import PhysiologyBuilder

def fetch_drug_from_pubchem(drug_name, csv_path='benchmark_drugs.csv'):
    try:
        import pubchempy as pcp
        compounds = pcp.get_compounds(drug_name, 'name')
        if not compounds: return None
        compound = compounds[0]
        smiles = getattr(compound, 'connectivity_smiles', None) or compound.smiles
        
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        fup = round(max(min(1 / (1 + 10**(0.45 * logp - 0.7)), 1.0), 0.01), 3)
        
        aromatic_atoms = sum([mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())])
        heavy_atoms = max(1, mol.GetNumHeavyAtoms())
        aromatic_proportion = aromatic_atoms / heavy_atoms
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable_bonds - 0.74 * aromatic_proportion
        solubility_mol_L = 10 ** logs
        cs_mg_L = round(solubility_mol_L * mw * 1000, 2)
        
        drug_dict = {
            'Drug': drug_name, 'MW_g_mol': mw, 'LogP': logp, 'pKa': 7.4, 
            'Type': 'Neutral', 'fup': fup, 'BP_Ratio': 0.60 if fup < 0.1 else 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Kp_Scalar': 1.0, 'kd_hr': 1.0, 
            'Cs_mg_L': cs_mg_L,
            'Vmax_hep_mg_hr': 50.0, 'Km_hep_mg_L': 5.0, 
            'Vmax_ren_mg_hr': 25.0, 'Km_ren_mg_L': 2.0, 
            'SMILES': smiles
        }
        
        try:
            df = pd.read_csv(csv_path)
        except FileNotFoundError:
            df = pd.DataFrame()
            
        if 'Drug' not in df.columns or drug_dict['Drug'] not in df['Drug'].values:
            df = pd.concat([df, pd.DataFrame([drug_dict])], ignore_index=True)
            df.to_csv(csv_path, index=False)

        return drug_dict
    except Exception as e:
        print(f"[!] PubChem lookup failed: {e}")
        return None

def load_metabolites_for_drug(drug_name, csv_path='metabolites.csv'):
    met_list = []
    try:
        df_mets = pd.read_csv(csv_path)
        if 'Parent_Drug' in df_mets.columns:
            filtered = df_mets[df_mets['Parent_Drug'].str.lower() == drug_name.lower()]
            for _, row in filtered.iterrows():
                m_params = {
                    'Drug': row['Metabolite_Name'], 'MW_g_mol': float(row['MW_g_mol']), 
                    'LogP': float(row['LogP']), 'pKa': 7.4, 'Type': 'Neutral', 
                    'fup': float(row['fup']), 'BP_Ratio': 0.90, 'Kp_Scalar': 1.0,
                    'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0 / len(filtered),
                    'Cs_mg_L': 10000.0, 
                    'Vmax_hep_mg_hr': 20.0, 'Km_hep_mg_L': 5.0, 
                    'Vmax_ren_mg_hr': 10.0, 'Km_ren_mg_L': 2.0
                }
                met_list.append(m_params)
    except FileNotFoundError:
        df_mets = pd.DataFrame()
        
    if not met_list:
        m_params_save = {
            'Parent_Drug': drug_name, 'Metabolite_Name': f"{drug_name}_Metabolite",
            'MW_g_mol': 200.0, 'LogP': 1.5, 'pKa': 7.4, 'Type': 'Neutral', 
            'fup': 0.15, 'BP_Ratio': 0.90, 'Kp_Scalar': 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0,
            'Cs_mg_L': 10000.0,
            'Vmax_hep_mg_hr': 20.0, 'Km_hep_mg_L': 5.0, 
            'Vmax_ren_mg_hr': 10.0, 'Km_ren_mg_L': 2.0, 'SMILES': 'N/A'
        }
        if 'Metabolite_Name' not in df_mets.columns or m_params_save['Metabolite_Name'] not in df_mets['Metabolite_Name'].values:
            df_mets = pd.concat([df_mets, pd.DataFrame([m_params_save])], ignore_index=True)
            df_mets.to_csv(csv_path, index=False)

        m_params_save['Drug'] = m_params_save['Metabolite_Name']
        met_list.append(m_params_save)
    return met_list

def process_dynamic_events(current_t, params, events_df):
    if events_df.empty: return params
    active_events = events_df[(events_df['Time_hr'] <= current_t) & (events_df['Time_hr'] > current_t - 1.0)]
    
    for _, event in active_events.iterrows():
        e_type = event['Event_Type']
        val = event['Value']
        
        if e_type == 'Food' and val == 'High_Fat_Meal':
            print(f"\n[EVENT @ {current_t}h]: High-Fat Meal. Gastric emptying delayed.")
            params['kt_stomach'] = 0.5  
            params['V_fluid_stomach'] = 0.50 
        elif e_type == 'Food' and val == 'Fasting':
            print(f"\n[EVENT @ {current_t}h]: Return to Fasting State.")
            params['kt_stomach'] = 2.0
            params['V_fluid_stomach'] = 0.25
        elif e_type == 'DDI_Inhibitor' and val == 'Started':
            print(f"\n[EVENT @ {current_t}h]: CYP Inhibitor Started. Hepatic CL drops 80%.")
            params['Dynamic_CL_Scalar'] = 0.20 
        elif e_type == 'DDI_Inhibitor' and val == 'Stopped':
            print(f"\n[EVENT @ {current_t}h]: CYP Inhibitor Stopped. Hepatic CL returning to baseline.")
            params['Dynamic_CL_Scalar'] = 1.00 
        elif e_type == 'DDI_Inducer' and val == 'Started':
            print(f"\n[EVENT @ {current_t}h]: CYP Inducer Started. Hepatic CL increases 200%.")
            params['Dynamic_CL_Scalar'] = 2.00 
        elif e_type == 'DDI_Inducer' and val == 'Stopped':
            print(f"\n[EVENT @ {current_t}h]: CYP Inducer Stopped. Hepatic CL returning to baseline.")
            params['Dynamic_CL_Scalar'] = 1.00 
        elif e_type == 'Smoking' and val == 'Started':
            print(f"\n[EVENT @ {current_t}h]: Smoking Induced. Hepatic CL increases 50%.")
            params['Dynamic_CL_Scalar'] = 1.50
        elif e_type == 'Smoking' and val == 'Stopped':
            print(f"\n[EVENT @ {current_t}h]: Smoking Cessation. Hepatic CL returning to baseline.")
            params['Dynamic_CL_Scalar'] = 1.00
            
    return params

def oral_absorption_wbpbpk_odes(t, A, tissues, params_parent, kp_dict_parent, met_params_list, kp_dicts_met):
    n_tissues = len(tissues)
    n_mets = len(met_params_list)
    n_gut_comps = 6 
    
    dA_dt = np.zeros_like(A)
    
    A_solid = A[0 : n_gut_comps]
    A_diss  = A[n_gut_comps : 2 * n_gut_comps]
    
    kt = [params_parent.get('kt_stomach', 2.0), 4.0, 0.5, 0.5, 0.22, 0.07]
    v_fluid = [params_parent.get('V_fluid_stomach', 0.25), 0.05, 0.15, 0.10, 0.05, 0.01] 
    
    kd = params_parent.get('kd_hr', 1.0)
    cs_limit = params_parent.get('Cs_mg_L', 1000.0)
    
    ka_base = 1.2 
    ka_reg = [0.0, ka_base * 3.0, ka_base * 2.0, ka_base * 1.0, ka_base * 0.1, ka_base * 0.05]

    diss_rates = np.zeros(6)
    for i in range(6):
        c_local = A_diss[i] / v_fluid[i]
        driving_force = max(0.0, (cs_limit - c_local) / cs_limit)
        diss_rates[i] = kd * A_solid[i] * driving_force

    dA_dt[0] = -kt[0]*A_solid[0] - diss_rates[0]
    for i in range(1, 6):
        dA_dt[i] = kt[i-1]*A_solid[i-1] - kt[i]*A_solid[i] - diss_rates[i]

    total_absorption_rate = 0.0
    dA_dt[6] = diss_rates[0] - kt[0]*A_diss[0] - ka_reg[0]*A_diss[0]
    total_absorption_rate += ka_reg[0]*A_diss[0]
    
    for i in range(1, 6):
        dA_dt[6+i] = diss_rates[i] + kt[i-1]*A_diss[i-1] - kt[i]*A_diss[i] - ka_reg[i]*A_diss[i]
        total_absorption_rate += ka_reg[i]*A_diss[i]

    parent_start_idx = 2 * n_gut_comps
    A_parent = A[parent_start_idx : parent_start_idx + n_tissues]
    
    parent_exc_idx = parent_start_idx + n_tissues
    idx_parent_feces = parent_exc_idx + 1
    idx_parent_abs_start = idx_parent_feces + 1
    idx_parent_biliary = idx_parent_abs_start + n_gut_comps
    idx_parent_renal_act = idx_parent_biliary + 1
    idx_parent_renal_pass = idx_parent_renal_act + 1
    idx_parent_metabolized = idx_parent_renal_pass + 1

    A_dict_p = {tissues[i]: A_parent[i] for i in range(n_tissues)}
    dA_dict_p = {name: 0.0 for name in tissues}
    
    C_art_pool_p = A_dict_p['arterial_blood'] / (params_parent['V_blood'] * 0.33)
    C_ven_out_p = {}
    parallel_tissues = [t for t in tissues if t not in ['venous_blood', 'arterial_blood', 'lungs']]
    
    for tissue in parallel_tissues:
        C_t = A_dict_p[tissue] / params_parent[f'V_{tissue}']
        C_ven_out_p[tissue] = C_t / (kp_dict_parent[tissue.capitalize()] / params_parent['BP_Ratio'])
            
    sum_parallel_Q = sum([params_parent[f'Q_{t}'] for t in parallel_tissues])
    Q_missing = params_parent['Q_lungs'] - sum_parallel_Q
    
    C_ven_pool_p = A_dict_p['venous_blood'] / (params_parent['V_blood'] * 0.67) 
    sum_ven_ret_p = sum([params_parent[f'Q_{t}'] * C_ven_out_p[t] for t in parallel_tissues])
    dA_dict_p['venous_blood'] = sum_ven_ret_p + (Q_missing * C_art_pool_p) - (params_parent['Q_lungs'] * C_ven_pool_p)
    
    C_lung_ven_p = (A_dict_p['lungs'] / params_parent['V_lungs']) / (kp_dict_parent['Lungs'] / params_parent['BP_Ratio'])
    dA_dict_p['lungs'] = params_parent['Q_lungs'] * (C_ven_pool_p - C_lung_ven_p)
    dA_dict_p['arterial_blood'] = params_parent['Q_lungs'] * (C_lung_ven_p - C_art_pool_p)
    
    total_metabolism_rate = 0.0
    renal_excretion_rate_p = 0.0
    
    for tissue in parallel_tissues:
        rate_in = params_parent[f'Q_{tissue}'] * C_art_pool_p
        rate_out = params_parent[f'Q_{tissue}'] * C_ven_out_p[tissue]
        
        if tissue == 'liver':
            mppgl = 45.0 
            liver_weight_g = params_parent['V_liver'] * 1000 
            C_unbound_liver = (A_dict_p['liver'] / params_parent['V_liver']) * (params_parent['fup'] / kp_dict_parent['Liver'])
            
            cl_int_total = (params_parent['Hepatic_CLint_uL_min_mg'] * mppgl * liver_weight_g * 60) / 1e6
            cl_scalar = params_parent.get('Hepatic_CL_Scalar', 1.0) * params_parent.get('Dynamic_CL_Scalar', 1.0)
            total_metabolism_rate = (cl_int_total * cl_scalar) * C_unbound_liver
            
            vmax_hep = params_parent.get('Vmax_hep_mg_hr', 0.0)
            km_hep = params_parent.get('Km_hep_mg_L', 1.0)
            active_biliary_efflux = (vmax_hep * C_unbound_liver) / (km_hep + C_unbound_liver)
            
            dA_dt[7] += active_biliary_efflux 
            dA_dt[idx_parent_biliary] = active_biliary_efflux 

            dA_dict_p[tissue] = rate_in - rate_out + total_absorption_rate - total_metabolism_rate - active_biliary_efflux
        elif tissue == 'kidney':
            C_plasma = C_art_pool_p / params_parent['BP_Ratio']
            renal_cl_scalar = params_parent.get('Renal_CL_Scalar', 1.0)
            passive_filtration = ((7.2 * params_parent['fup']) * C_plasma) * renal_cl_scalar
            
            vmax_ren = params_parent.get('Vmax_ren_mg_hr', 0.0)
            km_ren = params_parent.get('Km_ren_mg_L', 1.0)
            active_renal_secretion = (vmax_ren * C_plasma) / (km_ren + C_plasma)
            
            dA_dt[idx_parent_renal_act] = active_renal_secretion
            dA_dt[idx_parent_renal_pass] = passive_filtration
            
            renal_excretion_rate_p = passive_filtration + active_renal_secretion
            dA_dict_p[tissue] = rate_in - rate_out - renal_excretion_rate_p
        else:
            dA_dict_p[tissue] = rate_in - rate_out
            
    for i in range(n_tissues):
        dA_dt[parent_start_idx + i] = dA_dict_p[tissues[i]]

    dA_dt[parent_exc_idx] = renal_excretion_rate_p + total_metabolism_rate
    dA_dt[idx_parent_feces] = kt[5]*A_solid[5] + kt[5]*A_diss[5]
    dA_dt[idx_parent_metabolized] = total_metabolism_rate
    
    for i in range(n_gut_comps):
        dA_dt[idx_parent_abs_start + i] = ka_reg[i]*A_diss[i]

    # 3. METABOLITE TISSUE DISTRIBUTION
    for m_idx in range(n_mets):
        m_params = met_params_list[m_idx]
        kp_dict_m = kp_dicts_met[m_idx]
        
        m_start_idx = idx_parent_metabolized + 1 + m_idx * (n_tissues + 1)
        A_met = A[m_start_idx : m_start_idx + n_tissues]
        A_dict_m = {tissues[i]: A_met[i] for i in range(n_tissues)}
        dA_dict_m = {name: 0.0 for name in tissues}
        
        C_art_pool_m = A_dict_m['arterial_blood'] / (params_parent['V_blood'] * 0.33)
        C_ven_out_m = {}
        for tissue in parallel_tissues:
            C_t = A_dict_m[tissue] / params_parent[f'V_{tissue}']
            C_ven_out_m[tissue] = C_t / (kp_dict_m[tissue.capitalize()] / m_params['BP_Ratio'])
                
        C_ven_pool_m = A_dict_m['venous_blood'] / (params_parent['V_blood'] * 0.67) 
        sum_ven_ret_m = sum([params_parent[f'Q_{t}'] * C_ven_out_m[t] for t in parallel_tissues])
        dA_dict_m['venous_blood'] = sum_ven_ret_m + (Q_missing * C_art_pool_m) - (params_parent['Q_lungs'] * C_ven_pool_m)
        
        C_lung_ven_m = (A_dict_m['lungs'] / params_parent['V_lungs']) / (kp_dict_m['Lungs'] / m_params['BP_Ratio'])
        dA_dict_m['lungs'] = params_parent['Q_lungs'] * (C_ven_pool_m - C_lung_ven_m)
        dA_dict_m['arterial_blood'] = params_parent['Q_lungs'] * (C_lung_ven_m - C_art_pool_m)
        
        mw_ratio = m_params['MW_g_mol'] / params_parent['MW_g_mol']
        formed_this_metabolite = total_metabolism_rate * mw_ratio * m_params['Formation_Fraction']
        renal_excretion_rate_m = 0.0

        for tissue in parallel_tissues:
            rate_in = params_parent[f'Q_{tissue}'] * C_art_pool_m
            rate_out = params_parent[f'Q_{tissue}'] * C_ven_out_m[tissue]
            
            if tissue == 'liver':
                C_unbound_m_liver = (A_dict_m['liver'] / params_parent['V_liver']) * (m_params['fup'] / kp_dict_m['Liver'])
                vmax_hep_m = m_params.get('Vmax_hep_mg_hr', 0.0)
                km_hep_m = m_params.get('Km_hep_mg_L', 1.0)
                met_biliary_efflux = (vmax_hep_m * C_unbound_m_liver) / (km_hep_m + C_unbound_m_liver)
                dA_dict_m[tissue] = rate_in - rate_out + formed_this_metabolite - met_biliary_efflux
            elif tissue == 'kidney':
                C_plasma_m = C_art_pool_m / m_params['BP_Ratio']
                passive_m = 4.0 * m_params['fup'] * C_plasma_m
                vmax_ren_m = m_params.get('Vmax_ren_mg_hr', 0.0)
                km_ren_m = m_params.get('Km_ren_mg_L', 1.0)
                active_ren_m = (vmax_ren_m * C_plasma_m) / (km_ren_m + C_plasma_m)
                renal_excretion_rate_m = passive_m + active_ren_m
                dA_dict_m[tissue] = rate_in - rate_out - renal_excretion_rate_m
            else:
                dA_dict_m[tissue] = rate_in - rate_out
                
        for i in range(n_tissues):
            dA_dt[m_start_idx + i] = dA_dict_m[tissues[i]]
            
        met_exc_idx = m_start_idx + n_tissues
        dA_dt[met_exc_idx] = renal_excretion_rate_m
            
    return dA_dt

def run_simulation(drug_name, dose_mg, duration_days, interval_hours, weight_kg, cirrhosis, renal, events_csv, save_png=True, drug_csv='benchmark_drugs.csv', met_csv='metabolites.csv'):
    import matplotlib
    matplotlib.use('Agg') 
    import matplotlib.pyplot as plt

    try:
        df_drugs = pd.read_csv(drug_csv, index_col='Drug')
    except FileNotFoundError:
        df_drugs = pd.DataFrame()

    if drug_name in df_drugs.index:
        drug_params = df_drugs.loc[drug_name].to_dict()
        drug_params['Drug'] = drug_name
    else:
        drug_params = fetch_drug_from_pubchem(drug_name, drug_csv)
        if not drug_params: return None

    met_params_list = load_metabolites_for_drug(drug_name, met_csv)
    
    # Init Physiology Builder
    builder = PhysiologyBuilder()
    physio = builder.build_system(weight_kg=weight_kg, cirrhosis_status=cirrhosis, renal_status=renal)
    params = {**physio, **drug_params}
    
    params['kt_stomach'] = 2.0
    params['V_fluid_stomach'] = 0.25
    params['Dynamic_CL_Scalar'] = 1.0 
    
    kp_dict = calculate_tissue_kp(drug_params)
    kp_dicts_met = [calculate_tissue_kp(m) for m in met_params_list]
    
    tissues = ['venous_blood', 'lungs', 'arterial_blood', 'adipose', 'liver', 'muscle', 'kidney', 'skin', 'brain', 'heart', 'bone']
    n_tissues = len(tissues)
    n_mets = len(met_params_list)
    n_gut_comps = 6 
    
    total_compartments = (2 * n_gut_comps) + n_tissues + 1 + 1 + n_gut_comps + 4 + n_mets * (n_tissues + 1)
    
    dosing_hours = max(72.0, duration_days * 24.0) if duration_days > 0 else 72.0
    dosing_times = np.arange(0, duration_days * 24.0, interval_hours) if duration_days > 0 else [0.0]

    current_t = 0.0
    A_current = np.zeros(total_compartments)
    all_t, all_y = [], []

    # Load Dynamic Events
    try:
        events_df = pd.read_csv(events_csv)
    except FileNotFoundError:
        events_df = pd.DataFrame(columns=['Time_hr', 'Event_Type', 'Value'])

    log_filename = f"{drug_name}_PK.log"
    class Tee:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w")
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
            self.log.flush()
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    tee = Tee(log_filename)
    sys.stdout = tee

    print(f"\n[+] Initialization: PopPK Profile (Weight: {weight_kg}kg | Cirrhosis: {cirrhosis} | Renal: {renal})")
    
    # Dynamic Integration Loop
    while current_t < dosing_hours:
        step_duration = 1.0 
        if any(np.isclose(dosing_times, current_t, atol=1e-3)):
            A_current[0] += dose_mg
            
        params = process_dynamic_events(current_t, params, events_df)
        
        t_eval_seg = np.linspace(current_t, current_t + step_duration, 10)
        sol = solve_ivp(oral_absorption_wbpbpk_odes, (current_t, current_t + step_duration), A_current, 
                        t_eval=t_eval_seg, args=(tissues, params, kp_dict, met_params_list, kp_dicts_met), method='Radau')
        
        all_t.extend(sol.t[:-1]) 
        all_y.append(sol.y[:, :-1])
        A_current = sol.y[:, -1]
        current_t += step_duration

    temp_y = np.hstack(all_y)
    venous_idx = (2 * n_gut_comps) + tissues.index('venous_blood')
    
    plasma_conc = temp_y[venous_idx, :] / (params['V_blood'] * 0.67 * params['BP_Ratio'])
    c_max = np.max(plasma_conc)
    target_threshold = c_max / 100.0
    tail_duration = 1000.0
    t_end_tail = current_t + tail_duration
    
    def concentration_threshold_event(t, A, *args):
        venous_mass = A[venous_idx]
        conc = venous_mass / (params['V_blood'] * 0.67 * params['BP_Ratio'])
        return conc - target_threshold
    concentration_threshold_event.terminal = True
    concentration_threshold_event.direction = -1

    sol_tail = solve_ivp(oral_absorption_wbpbpk_odes, (current_t, t_end_tail), A_current,
                         t_eval=np.linspace(current_t, t_end_tail, 200),
                         args=(tissues, params, kp_dict, met_params_list, kp_dicts_met),
                         method='Radau', events=concentration_threshold_event)

    all_t.extend(sol_tail.t)
    all_y.append(sol_tail.y)

    sol_t = np.array(all_t)
    sol_y = np.hstack(all_y)

    # -------------------------------------------------------------
    # NON-COMPARTMENTAL ANALYSIS (NCA) PK PARAMETERS
    # -------------------------------------------------------------
    
    idx_feces = (2 * n_gut_comps) + n_tissues + 1
    idx_abs_start = idx_feces + 1
    idx_biliary = idx_abs_start + 6
    idx_renal_act = idx_biliary + 1
    idx_renal_pass = idx_renal_act + 1
    idx_metabolized = idx_renal_pass + 1

    final_feces = sol_y[idx_feces, -1]
    final_abs = [sol_y[idx_abs_start + i, -1] for i in range(n_gut_comps)]
    total_abs = sum(final_abs)
    
    tot_biliary = sol_y[idx_biliary, -1]
    tot_ren_act = sol_y[idx_renal_act, -1]
    tot_ren_pass = sol_y[idx_renal_pass, -1]
    tot_metabolized = sol_y[idx_metabolized, -1]

    total_dose_administered = dose_mg * len(dosing_times)

    full_plasma_conc = sol_y[venous_idx, :] / (params['V_blood'] * 0.67 * params['BP_Ratio'])
    auc_0_t = np.trapz(full_plasma_conc, sol_t)

    valid_idx = full_plasma_conc > 1e-6 
    if np.sum(valid_idx) > 10:
        tail_t = sol_t[valid_idx][-50:]
        tail_c = full_plasma_conc[valid_idx][-50:]
        try:
            slope, _ = np.polyfit(tail_t, np.log(tail_c), 1)
            lambda_z = -slope
        except:
            lambda_z = 0
    else:
        lambda_z = 0

    t_half = np.log(2) / lambda_z if lambda_z > 0 else float('inf')
    auc_inf = auc_0_t + (full_plasma_conc[-1] / lambda_z if lambda_z > 0 else 0)

    systemic_elimination = tot_metabolized + tot_ren_act + tot_ren_pass
    F_bioavailability = systemic_elimination / total_dose_administered if total_dose_administered > 0 else 0
    
    cl_total = systemic_elimination / auc_inf if auc_inf > 0 else 0
    cl_hepatic = tot_metabolized / auc_inf if auc_inf > 0 else 0
    cl_renal = (tot_ren_act + tot_ren_pass) / auc_inf if auc_inf > 0 else 0
    
    cl_apparent = total_dose_administered / auc_inf if auc_inf > 0 else 0

    c_max_actual = np.max(full_plasma_conc)
    t_max_actual = sol_t[np.argmax(full_plasma_conc)]

    print(f"\n[*] Target: 1/100th Cmax = {target_threshold:.4f} mg/L | Total Time: {sol_t[-1]:.1f}h")
    print("\n" + "="*40)
    print(f"  PK PREDICTION RESULTS: {drug_name.upper()}")
    print("="*40)
    
    print("\n  NON-COMPARTMENTAL PK PARAMETERS (NCA)")
    print("="*40)
    print(f"  Cmax:                  {c_max_actual:>10.4f} mg/L")
    print(f"  Tmax:                  {t_max_actual:>10.2f} h")
    print(f"  AUC (0-inf):           {auc_inf:>10.2f} mg*h/L")
    print(f"  t1/2 (Terminal):       {t_half:>10.2f} h")
    print(f"  Fraction Absorbed (Fa):{(total_abs/total_dose_administered)*100:>10.2f} %")
    print(f"  Bioavailability (F):   {F_bioavailability * 100:>10.2f} %")
    print(f"  Apparent F:            {F_bioavailability * 100:>10.2f} %")
    print(f"  CL_hepatic:            {cl_hepatic:>10.4f} L/h")
    print(f"  CL_renal:              {cl_renal:>10.4f} L/h")
    print(f"  CL_total (Systemic):   {cl_total:>10.4f} L/h")
    print(f"  Apparent CL (CL/F):    {cl_apparent:>10.4f} L/h")
    
    print("\n" + "="*40)
    print("  GIT ABSORPTION & EXCRETION SUMMARY")
    print("="*40)
    print(f"Total Amount Absorbed: {total_abs:.4f} mg (Includes EHC re-absorption)")
    print(f"Amount Excreted (Feces): {final_feces:.4f} mg")
    print("-" * 40)
    print(f"  Stomach:  {final_abs[0]:>10.4f} mg")
    print(f"  Duodenum: {final_abs[1]:>10.4f} mg")
    print(f"  Jejunum:  {final_abs[2]:>10.4f} mg")
    print(f"  Ileum:    {final_abs[3]:>10.4f} mg")
    print(f"  Cecum:    {final_abs[4]:>10.4f} mg")
    print(f"  Colon:    {final_abs[5]:>10.4f} mg")
    
    print("\n" + "="*40)
    print("  ACTIVE TRANSPORT & METABOLISM FLUX")
    print("="*40)
    print(f"  Hepatic Mass Metabolized:   {tot_metabolized:>10.4f} mg")
    print(f"  Biliary Efflux (EHC route): {tot_biliary:>10.4f} mg")
    print(f"  Renal Active Secretion:     {tot_ren_act:>10.4f} mg")
    print(f"  Renal Passive GFR:          {tot_ren_pass:>10.4f} mg")

    m_col_widths = []
    for m in met_params_list:
        m_name = m['Drug']
        col_w = max(16, len(m_name) + 6) 
        m_col_widths.append(col_w)

    header_str = f"  {'Tissue':<15} | {'Parent Max Mass':<16} | {'Parent Peak Conc':<18}"
    for idx, m in enumerate(met_params_list):
        m_name = m['Drug']
        col_w = m_col_widths[idx]
        header_str += f" | {m_name + '_Mass':<{col_w}} | {m_name + '_Conc':<{col_w}}"
    
    print("\n" + "="*len(header_str))
    print("  TISSUE ACCUMULATION SUMMARY")
    print("="*len(header_str))
    print(header_str)
    print("-" * len(header_str))

    for idx, tissue in enumerate(tissues):
        p_row_idx = (2 * n_gut_comps) + idx
        p_masses = sol_y[p_row_idx, :]
        p_max_mass = np.max(p_masses)
        vol = params.get(f'V_{tissue}', params['V_blood'] * 0.5)
        p_peak_conc = p_max_mass / vol
        row_output = f"  {tissue.capitalize():<15} | {p_max_mass:<16.4f} | {p_peak_conc:<18.4f}"

        for m_idx in range(n_mets):
            m_row_idx = idx_metabolized + 1 + m_idx * (n_tissues + 1) + idx
            m_masses = sol_y[m_row_idx, :]
            m_max_mass = np.max(m_masses)
            m_peak_conc = m_max_mass / vol
            col_w = m_col_widths[m_idx]
            row_output += f" | {m_max_mass:<{col_w}.4f} | {m_peak_conc:<{col_w}.4f}"
        print(row_output)
    print("=" * len(header_str))

    # -------------------------------------------------------------
    # CSV EXPORT: TIME VS CONCENTRATION
    # -------------------------------------------------------------
    csv_data = {'Time_hr': sol_t}
    
    for idx, tissue in enumerate(tissues):
        p_row_idx = (2 * n_gut_comps) + idx
        vol = params.get(f'V_{tissue}', params['V_blood'] * 0.5)
        if tissue == 'venous_blood':
            csv_data['Parent_Venous_Plasma_mg_L'] = full_plasma_conc
        elif tissue == 'arterial_blood':
            csv_data['Parent_Arterial_Plasma_mg_L'] = sol_y[p_row_idx, :] / (params['V_blood'] * 0.33 * params['BP_Ratio'])
        else:
            csv_data[f'Parent_{tissue.capitalize()}_mg_L'] = sol_y[p_row_idx, :] / vol

    for m_idx, m_params in enumerate(met_params_list):
        m_name = m_params['Drug']
        for idx, tissue in enumerate(tissues):
            m_row_idx = idx_metabolized + 1 + m_idx * (n_tissues + 1) + idx
            vol = params.get(f'V_{tissue}', params['V_blood'] * 0.5)
            if tissue == 'venous_blood':
                csv_data[f'{m_name}_Venous_Plasma_mg_L'] = sol_y[m_row_idx, :] / (params['V_blood'] * 0.67 * m_params['BP_Ratio'])
            elif tissue == 'arterial_blood':
                csv_data[f'{m_name}_Arterial_Plasma_mg_L'] = sol_y[m_row_idx, :] / (params['V_blood'] * 0.33 * m_params['BP_Ratio'])
            else:
                csv_data[f'{m_name}_{tissue.capitalize()}_mg_L'] = sol_y[m_row_idx, :] / vol

    df_csv = pd.DataFrame(csv_data)
    csv_filename = f"{drug_name}_Time_Course.csv"
    df_csv.to_csv(csv_filename, index=False)
    print(f"\n[+] Time-course concentration data exported successfully to {csv_filename}")

    if save_png:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(sol_t, np.maximum(full_plasma_conc, 1e-6), label=f'Venous Plasma (Parent)', color='#0088cc', lw=2)
        
        met_colors = ['#ff6600', '#009966', '#cc00cc', '#d62728'] 
        for m_idx, m_params in enumerate(met_params_list):
            m_venous_idx = idx_metabolized + 1 + m_idx * (n_tissues + 1) + tissues.index('venous_blood')
            m_plasma_conc = sol_y[m_venous_idx, :] / (params['V_blood'] * 0.67 * m_params['BP_Ratio'])
            c = met_colors[m_idx % len(met_colors)]
            ax.plot(sol_t, np.maximum(m_plasma_conc, 1e-6), label=f"Venous Plasma ({m_params['Drug']})", color=c, lw=2, linestyle='-')

        ax.axhline(target_threshold, color='red', linestyle='--', label='1/100th Parent Cmax')
        
        ax.set_yscale('log')
        ax.set_ylim(bottom=min(1e-4, target_threshold / 5.0))
        ax.set_title(f'PK Profile (Phase 4 QSP): {drug_name}')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (mg/L)')
        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        ax.grid(True, alpha=0.3, which="both", ls="--")
        
        png_name = f"{drug_name}_Phase4_PopPK_Profile.png"
        plt.savefig(png_name, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"[+] Graph saved successfully as {png_name}")

    sys.stdout = sys.stdout.terminal
    return sol_t, sol_y

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master Controller PBPK Simulation (Phase 4)")
    parser.add_argument('--drug', type=str, default='Efavirenz')
    parser.add_argument('--dose', type=float, default=100.0)
    parser.add_argument('--duration_days', type=float, default=30.0)
    parser.add_argument('--interval_hours', type=float, default=24.0)
    parser.add_argument('--weight', type=float, default=70.0)
    parser.add_argument('--cirrhosis', type=str, default='None')
    parser.add_argument('--renal', type=str, default='Healthy')
    parser.add_argument('--events_csv', type=str, default='patient_events.csv')
    parser.add_argument('--save_png', action='store_true')
    parser.add_argument('--drug_csv', type=str, default='benchmark_drugs.csv')
    parser.add_argument('--met_csv', type=str, default='metabolites.csv')
    args = parser.parse_args()
    
    run_simulation(args.drug, args.dose, args.duration_days, args.interval_hours, 
                   args.weight, args.cirrhosis, args.renal, args.events_csv,
                   args.save_png, args.drug_csv, args.met_csv)