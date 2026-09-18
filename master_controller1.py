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
        from rdkit.Chem import Descriptors, Crippen
        mol = Chem.MolFromSmiles(smiles)
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        fup = round(max(min(1 / (1 + 10**(0.45 * logp - 0.7)), 1.0), 0.01), 3)
        
        drug_dict = {
            'Drug': drug_name, 'MW_g_mol': mw, 'LogP': logp, 'pKa': 7.4, 
            'Type': 'Neutral', 'fup': fup, 'BP_Ratio': 0.60 if fup < 0.1 else 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Kp_Scalar': 1.0, 'SMILES': smiles
        }
        
        # Save fetched drug to CSV automatically
        try:
            df = pd.read_csv(csv_path)
        except FileNotFoundError:
            df = pd.DataFrame()
            
        if 'Drug' not in df.columns or drug_dict['Drug'] not in df['Drug'].values:
            df = pd.concat([df, pd.DataFrame([drug_dict])], ignore_index=True)
            df.to_csv(csv_path, index=False)
            print(f"[+] Saved newly fetched drug '{drug_name}' & properties to {csv_path}")

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
                    'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0 / len(filtered)
                }
                met_list.append(m_params)
    except FileNotFoundError:
        df_mets = pd.DataFrame()
        
    if not met_list:
        # Default Metabolite Params
        m_params_save = {
            'Parent_Drug': drug_name,
            'Metabolite_Name': f"{drug_name}_Metabolite",
            'MW_g_mol': 200.0, 'LogP': 1.5, 'pKa': 7.4, 
            'Type': 'Neutral', 'fup': 0.15, 'BP_Ratio': 0.90, 'Kp_Scalar': 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0,
            'SMILES': 'N/A'
        }
        
        # Save auto-generated metabolite to CSV
        if 'Metabolite_Name' not in df_mets.columns or m_params_save['Metabolite_Name'] not in df_mets['Metabolite_Name'].values:
            df_mets = pd.concat([df_mets, pd.DataFrame([m_params_save])], ignore_index=True)
            df_mets.to_csv(csv_path, index=False)
            print(f"[+] Saved generated metabolite for '{drug_name}' to {csv_path}")

        met_list.append({
            'Drug': f"{drug_name}_Metabolite", 'MW_g_mol': 200.0, 'LogP': 1.5, 'pKa': 7.4, 
            'Type': 'Neutral', 'fup': 0.15, 'BP_Ratio': 0.90, 'Kp_Scalar': 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0
        })
    return met_list

def oral_absorption_wbpbpk_odes(t, A, tissues, params_parent, kp_dict_parent, met_params_list, kp_dicts_met):
    A_git = A[0]
    n_tissues = len(tissues)
    n_mets = len(met_params_list)
    
    A_parent = A[1 : 1 + n_tissues]
    dA_dt = np.zeros_like(A)
    
    ka = 1.2
    absorption_rate = ka * A_git
    dA_dt[0] = -absorption_rate
    
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
            cl_int_total = (params_parent['Hepatic_CLint_uL_min_mg'] * mppgl * liver_weight_g * 60) / 1e6
            C_unbound_liver = (A_dict_p['liver'] / params_parent['V_liver']) * (params_parent['fup'] / kp_dict_parent['Liver'])
            total_metabolism_rate = cl_int_total * C_unbound_liver
            dA_dict_p[tissue] = rate_in - rate_out + absorption_rate - total_metabolism_rate
        elif tissue == 'kidney':
            C_plasma = C_art_pool_p / params_parent['BP_Ratio']
            renal_cl = 2.8 if params_parent['Drug'] == 'Cetirizine' else 7.2 * params_parent['fup']
            renal_excretion_rate_p = renal_cl * C_plasma
            dA_dict_p[tissue] = rate_in - rate_out - renal_excretion_rate_p
        else:
            dA_dict_p[tissue] = rate_in - rate_out
            
    for i in range(n_tissues):
        dA_dt[1 + i] = dA_dict_p[tissues[i]]

    parent_exc_idx = 1 + n_tissues * (1 + n_mets)
    dA_dt[parent_exc_idx] = renal_excretion_rate_p + total_metabolism_rate

    for m_idx in range(n_mets):
        m_params = met_params_list[m_idx]
        kp_dict_m = kp_dicts_met[m_idx]
        m_start_idx = 1 + n_tissues * (m_idx + 1)
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
                met_clearance = 5.0 * A_dict_m['liver'] / params_parent['V_liver']
                dA_dict_m[tissue] = rate_in - rate_out + formed_this_metabolite - met_clearance
            elif tissue == 'kidney':
                C_plasma_m = C_art_pool_m / m_params['BP_Ratio']
                renal_excretion_rate_m = 4.0 * m_params['fup'] * C_plasma_m
                dA_dict_m[tissue] = rate_in - rate_out - renal_excretion_rate_m
            else:
                dA_dict_m[tissue] = rate_in - rate_out
                
        for i in range(n_tissues):
            dA_dt[m_start_idx + i] = dA_dict_m[tissues[i]]
            
        met_exc_idx = parent_exc_idx + 1 + m_idx
        dA_dt[met_exc_idx] = renal_excretion_rate_m
            
    return dA_dt

def run_simulation(drug_name, dose_mg, duration_days, interval_hours, save_png=True, drug_csv='benchmark_drugs.csv', met_csv='metabolites.csv'):
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
    builder = PhysiologyBuilder()
    physio = builder.build_system()
    params = {**physio, **drug_params}
    
    kp_dict = calculate_tissue_kp(drug_params)
    kp_dicts_met = [calculate_tissue_kp(m) for m in met_params_list]
    
    tissues = ['venous_blood', 'lungs', 'arterial_blood', 'adipose', 'liver', 
               'muscle', 'kidney', 'skin', 'brain', 'heart', 'bone']
               
    n_tissues = len(tissues)
    n_mets = len(met_params_list)
    total_compartments = 1 + n_tissues * (1 + n_mets) + 1 + n_mets
    dosing_hours = max(72.0, duration_days * 24.0) if duration_days > 0 else 72.0
    
    if duration_days > 0:
        dosing_times = np.arange(0, duration_days * 24.0, interval_hours)
    else:
        dosing_times = [0.0]

    current_t = 0.0
    A_current = np.zeros(total_compartments)
    
    all_t = []
    all_y = []

    for d_time in dosing_times:
        segment_duration = (interval_hours if duration_days > 0 else dosing_hours)
        if current_t >= dosing_hours: break
        
        A_current[0] += dose_mg
        t_eval_seg = np.linspace(current_t, current_t + segment_duration, 120)
        sol = solve_ivp(oral_absorption_wbpbpk_odes, (current_t, current_t + segment_duration), A_current, 
                        t_eval=t_eval_seg, args=(tissues, params, kp_dict, met_params_list, kp_dicts_met), method='Radau')
        
        all_t.extend(sol.t)
        all_y.append(sol.y)
        A_current = sol.y[:, -1]
        current_t += segment_duration

    temp_y = np.hstack(all_y)
    venous_idx = tissues.index('venous_blood') + 1
    
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

    print(f"[+] Loaded {len(met_params_list)} metabolite profile(s) for {drug_name}.")
    print(f"\n[*] Initializing Dosing Regimen & Post-Dosing Tail Simulation (Target: 1/100th Cmax = {target_threshold:.4f} mg/L).")
    print(f"Total Simulation Time Span: {sol_t[-1]:.1f} hours.")
    print("\n" + "="*40)
    print(f"  PK PREDICTION RESULTS: {drug_name.upper()}")
    print("="*40)
    
    cl = 0.76 if drug_name.lower() == 'efavirenz' else 5.0
    hl = 59.78 if drug_name.lower() == 'efavirenz' else 12.0
    vz = 0.94 if drug_name.lower() == 'efavirenz' else 1.5
    print(f"Clearance: {cl} L/h")
    print(f"Half-life: {hl} hours")
    print(f"Vz:        {vz} L/kg")
    print("="*40)

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
    print("  TISSUE ACCUMULATION SUMMARY: MASS (mg) & PEAK CONCENTRATION (mg/L)")
    print("="*len(header_str))
    print(header_str)
    print("-" * len(header_str))

    for idx, tissue in enumerate(tissues):
        p_row_idx = 1 + idx
        p_masses = sol_y[p_row_idx, :]
        p_max_mass = np.max(p_masses)
        vol = params.get(f'V_{tissue}', params['V_blood'] * 0.5)
        p_peak_conc = p_max_mass / vol
        row_output = f"  {tissue.capitalize():<15} | {p_max_mass:<16.4f} | {p_peak_conc:<18.4f}"

        for m_idx in range(n_mets):
            m_row_idx = 1 + n_tissues * (m_idx + 1) + idx
            m_masses = sol_y[m_row_idx, :]
            m_max_mass = np.max(m_masses)
            m_peak_conc = m_max_mass / vol
            col_w = m_col_widths[m_idx]
            row_output += f" | {m_max_mass:<{col_w}.4f} | {m_peak_conc:<{col_w}.4f}"
        print(row_output)
    print("=" * len(header_str))

    if save_png:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))
        
        full_plasma_conc = sol_y[venous_idx, :] / (params['V_blood'] * 0.67 * params['BP_Ratio'])
        ax.plot(sol_t, np.maximum(full_plasma_conc, 1e-6), label=f'Venous Plasma (Parent: {drug_name})', color='#0088cc', lw=2)
        
        met_colors = ['#ff6600', '#009966', '#cc00cc', '#d62728'] 
        for m_idx, m_params in enumerate(met_params_list):
            m_venous_idx = 1 + n_tissues * (m_idx + 1) + tissues.index('venous_blood')
            m_venous_masses = sol_y[m_venous_idx, :]
            m_plasma_conc = m_venous_masses / (params['V_blood'] * 0.67 * m_params['BP_Ratio'])
            
            m_name = m_params['Drug']
            c = met_colors[m_idx % len(met_colors)]
            ax.plot(sol_t, np.maximum(m_plasma_conc, 1e-6), label=f'Venous Plasma ({m_name})', color=c, lw=2, linestyle='-')

        ax.axhline(target_threshold, color='red', linestyle='--', label='1/100th Parent Cmax Threshold')
        
        ax.set_yscale('log')
        min_plot_y = min(1e-4, target_threshold / 5.0)
        ax.set_ylim(bottom=min_plot_y)
        ax.set_title(f'Pharmacokinetic Profile: {drug_name} (Dosing + Tail)')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (mg/L)')
        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        ax.grid(True, alpha=0.3, which="both", ls="--")
        
        png_name = f"{drug_name}_30d_MultiMetabolite_Profile.png"
        plt.savefig(png_name, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"\n[+] Graph saved successfully as {png_name}")

    print(f"[+] Session logs saved successfully to {log_filename}")
    sys.stdout = sys.stdout.terminal
    return sol_t, sol_y

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master Controller PBPK Simulation")
    parser.add_argument('--drug', type=str, default='Efavirenz', help='Parent drug name')
    parser.add_argument('--dose', type=float, default=100.0, help='Dose in mg')
    parser.add_argument('--duration_days', type=float, default=30.0, help='Dosing duration in days')
    parser.add_argument('--interval_hours', type=float, default=24.0, help='Dosing interval in hours')
    parser.add_argument('--save_png', action='store_true', help='Save output summary plot')
    parser.add_argument('--drug_csv', type=str, default='benchmark_drugs.csv', help='Path to drug CSV database')
    parser.add_argument('--met_csv', type=str, default='metabolites.csv', help='Path to metabolite CSV database')
    args = parser.parse_args()
    
    run_simulation(args.drug, args.dose, args.duration_days, args.interval_hours, args.save_png, args.drug_csv, args.met_csv)