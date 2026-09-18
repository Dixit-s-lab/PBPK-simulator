import matplotlib
matplotlib.use('TkAgg', force=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
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
        m_params_save = {
            'Parent_Drug': drug_name,
            'Metabolite_Name': f"{drug_name}_Metabolite",
            'MW_g_mol': 200.0, 'LogP': 1.5, 'pKa': 7.4, 
            'Type': 'Neutral', 'fup': 0.15, 'BP_Ratio': 0.90, 'Kp_Scalar': 1.0,
            'Hepatic_CLint_uL_min_mg': 1.0, 'Formation_Fraction': 1.0,
            'SMILES': 'N/A'
        }
        
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

def run_animator(drug_name, dose_mg, duration_days, interval_hours, anim_type='body', save_gif=False, drug_csv='benchmark_drugs.csv', met_csv='metabolites.csv'):
    try:
        df_drugs = pd.read_csv(drug_csv, index_col='Drug')
    except FileNotFoundError:
        df_drugs = pd.DataFrame()

    if drug_name in df_drugs.index:
        drug_params = df_drugs.loc[drug_name].to_dict()
        drug_params['Drug'] = drug_name
    else:
        drug_params = fetch_drug_from_pubchem(drug_name, drug_csv)
        if not drug_params: return

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
        t_eval_seg = np.linspace(current_t, current_t + segment_duration, 100)
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

    def concentration_threshold_event(t, A, *args):
        venous_mass = A[venous_idx]
        conc = venous_mass / (params['V_blood'] * 0.67 * params['BP_Ratio'])
        return conc - target_threshold
    concentration_threshold_event.terminal = True
    concentration_threshold_event.direction = -1

    sol_tail = solve_ivp(oral_absorption_wbpbpk_odes, (current_t, current_t + 1000.0), A_current,
                         t_eval=np.linspace(current_t, current_t + 1000.0, 150),
                         args=(tissues, params, kp_dict, met_params_list, kp_dicts_met),
                         method='Radau', events=concentration_threshold_event)

    all_t.extend(sol_tail.t)
    all_y.append(sol_tail.y)

    sol_t = np.array(all_t)
    sol_y = np.hstack(all_y)

    animations = []

    if anim_type in ['body', 'both']:
        fig_body, ax_body = plt.subplots(figsize=(7, 9))
        fig_body.patch.set_facecolor('white')
        ax_body.set_facecolor('white')
        plt.subplots_adjust(top=0.84, bottom=0.05, left=0.05, right=0.95)

        head = patches.Ellipse((0, 8.2), 1.8, 2.2, facecolor='none', edgecolor='black', linewidth=1.5)
        torso = patches.Rectangle((-1.3, 1.5), 2.6, 5.2, facecolor='none', edgecolor='black', linewidth=1.5)
        l_arm = patches.Rectangle((-2.4, 2.0), 0.8, 4.0, facecolor='none', edgecolor='black', linewidth=1.2)
        r_arm = patches.Rectangle((1.6, 2.0), 0.8, 4.0, facecolor='none', edgecolor='black', linewidth=1.2)
        l_leg = patches.Rectangle((-1.1, -3.2), 0.9, 4.8, facecolor='none', edgecolor='black', linewidth=1.2)
        r_leg = patches.Rectangle((0.2, -3.2), 0.9, 4.8, facecolor='none', edgecolor='black', linewidth=1.2)

        for p in [head, torso, l_arm, r_arm, l_leg, r_leg]:
            ax_body.add_patch(p)

        organs_pos = {
            'mouth': (0, 7.3), 'git': (0.0, 5.15),
            'brain': (0, 8.2), 'lungs': (0, 6.2), 'heart': (0.4, 5.2),
            'liver': (-0.8, 4.0), 'kidney': (0.6, 3.1), 'adipose': (-0.8, 5.5),
            'muscle_left': (-2.0, 4.0), 'muscle_right': (2.0, 4.0),
            'bone_left': (-0.65, -0.8), 'bone_right': (0.65, -0.8),
            'excretion': (0.0, 2.3)
        }
        heart_x, heart_y = organs_pos['heart']

        for name, pos in organs_pos.items():
            if name not in ['heart', 'lungs', 'mouth', 'git', 'muscle_left', 'muscle_right', 'bone_left', 'bone_right', 'excretion']:
                ax_body.plot([heart_x, pos[0]], [heart_y, pos[1]], color='#cc0000', linewidth=2, alpha=0.4, solid_capstyle='round')
                ax_body.plot([pos[0], heart_x], [pos[1], heart_y], color='#0066cc', linewidth=2, alpha=0.4, solid_capstyle='round')

        ax_body.plot([heart_x, 0.0], [heart_y, 6.2], color='#0066cc', linewidth=2.5, alpha=0.5, linestyle='--')
        ax_body.plot([0.0, heart_x], [6.2, heart_y], color='#cc0000', linewidth=2.5, alpha=0.5, linestyle='--')
        ax_body.plot([-0.5, -0.8], [2.3, 4.0], color='#008000', linewidth=2, alpha=0.3, solid_capstyle='round')
        ax_body.plot([0.6, 0.0], [3.1, 2.3], color='#cc0000', linewidth=1.5, alpha=0.3, linestyle='--')

        git_oval = patches.Ellipse((0.0, 5.15), 0.6, 3.3, facecolor='#e2f0d9', edgecolor='black', linewidth=1.2, alpha=0.35)
        ax_body.add_patch(git_oval)

        git_bottom_x, git_bottom_y = 0.0, 5.15 - (3.3 / 2.0)
        liver_x, liver_y = -0.8, 4.0
        ax_body.plot([git_bottom_x, liver_x], [git_bottom_y, liver_y], color='#008000', linewidth=3, alpha=0.6, solid_capstyle='round')
        ax_body.text((git_bottom_x + liver_x) / 2 - 0.15, (git_bottom_y + liver_y) / 2 - 0.3, 'Portal Vein', color='black', fontsize=6, fontweight='bold', ha='center', va='center')

        organ_coords = {
            'mouth': (0, 7.3, 0.5, 0.3), 'brain': (0, 8.2, 1.4, 1.1),
            'lungs': (0, 6.2, 1.3, 0.7), 'heart': (0.4, 5.2, 0.6, 0.6),
            'liver': (-0.8, 4.0, 1.0, 0.6), 'kidney': (0.6, 3.1, 0.5, 0.5),
            'adipose': (-0.8, 5.5, 0.5, 0.7), 'muscle': (0, 4.0, 3.2, 1.5),
            'skin': (0, 4.2, 2.6, 5.2), 'bone': (0, -0.8, 1.5, 2.5),
            'excretion': (0.0, 2.3, 0.5, 0.4), 'venous_blood': (-0.3, 5.8, 0.3, 0.3),
            'arterial_blood': (0.3, 5.8, 0.3, 0.3)
        }

        patches_dict = {'git': git_oval}
        ax_body.text(0.0, 3.7, 'GIT', color='black', fontsize=7, fontweight='bold', ha='center', va='center')

        for tissue, (x, y, w, h) in organ_coords.items():
            if tissue in ['venous_blood', 'arterial_blood']: continue
            is_exc = (tissue == 'excretion')
            patch = patches.Ellipse((x, y), w, h, facecolor='none', edgecolor='black', linewidth=1.2)
            ax_body.add_patch(patch)
            patches_dict[tissue] = patch
            if tissue == 'muscle':
                ax_body.text(-2.0, 4.0, 'Muscle', color='black', fontsize=7, fontweight='bold', ha='center', va='center')
                ax_body.text(2.0, 4.0, 'Muscle', color='black', fontsize=7, fontweight='bold', ha='center', va='center')
            elif tissue == 'bone':
                ax_body.text(-0.65, -0.8, 'Bone', color='black', fontsize=7, fontweight='bold', ha='center', va='center')
                ax_body.text(0.65, -0.8, 'Bone', color='black', fontsize=7, fontweight='bold', ha='center', va='center')
            elif is_exc:
                ax_body.text(x, y, 'Bladder', color='black', fontsize=6, fontweight='bold', ha='center', va='center')
            else:
                ax_body.text(x, y, tissue.capitalize(), color='black', fontsize=7, fontweight='bold', ha='center', va='center')

        ax_body.text(1.1, 4.2, 'Skin', color='black', fontsize=8, fontweight='bold', ha='left', va='center')

        compounds_info = [{'name': f"Parent: {drug_name}", 'color': '#0088cc', 'start_idx': 1, 'exc_idx': 1 + n_tissues * (1 + n_mets)}]
        met_colors = ['#cc00cc', '#ff6600', '#009966', '#d62728']
        for idx, m_params in enumerate(met_params_list):
            c = met_colors[idx % len(met_colors)]
            compounds_info.append({
                'name': f"Metabolite: {m_params['Drug']}", 
                'color': c, 'start_idx': 1 + n_tissues * (idx + 1),
                'exc_idx': (1 + n_tissues * (1 + n_mets)) + 1 + idx
            })

        scatter_handles = []
        for comp in compounds_info:
            sc = ax_body.scatter([], [], s=16, c=comp['color'], alpha=0.45, edgecolors='none', label=comp['name'])
            scatter_handles.append(sc)

        ax_body.legend(loc='upper right', facecolor='white', edgecolor='black', labelcolor='black', fontsize=8)
        ax_body.set_xlim(-3.2, 3.2)
        ax_body.set_ylim(-4.5, 10.5)
        ax_body.axis('off')

        def update_body(frame):
            git_mass = max(0.0, sol_y[0, frame])
            total_body_mass = np.sum(sol_y[:, frame]) + 1e-6
            np.random.seed(42 + frame)
            total_particles_per_comp = 200

            for c_idx, comp in enumerate(compounds_info):
                start_idx, exc_idx = comp['start_idx'], comp['exc_idx']
                comp_tissue_masses = sol_y[start_idx : start_idx + n_tissues, frame]
                exc_mass = max(0.0, sol_y[exc_idx, frame])
                comp_total_mass = np.sum(comp_tissue_masses) + (git_mass if c_idx == 0 else 0) + exc_mass + 1e-6

                particle_x, particle_y = [], []
                for idx, tissue in enumerate(tissues):
                    m_val = max(0.0, comp_tissue_masses[idx])
                    fraction = m_val / comp_total_mass
                    n_pts = int(fraction * total_particles_per_comp)
                    if n_pts <= 0: continue

                    if tissue == 'muscle':
                        half_pts = n_pts // 2
                        px = np.concatenate([np.random.uniform(-2.4, -1.6, half_pts), np.random.uniform(1.6, 2.4, n_pts - half_pts)])
                        py = np.concatenate([np.random.uniform(2.0, 6.0, half_pts), np.random.uniform(2.0, 6.0, n_pts - half_pts)])
                    elif tissue == 'bone':
                        half_pts = n_pts // 2
                        px = np.concatenate([np.random.uniform(-0.75, -0.55, half_pts), np.random.uniform(0.55, 0.75, n_pts - half_pts)])
                        py = np.concatenate([np.random.uniform(-1.6, 0.0, half_pts), np.random.uniform(-1.6, 0.0, n_pts - half_pts)])
                    else:
                        coords = organ_coords.get(tissue, (0, 0, 1, 1))
                        px = np.random.uniform(coords[0] - coords[2]/2, coords[0] + coords[2]/2, n_pts)
                        py = np.random.uniform(coords[1] - coords[3]/2, coords[1] + coords[3]/2, n_pts)

                    particle_x.extend(px); particle_y.extend(py)

                if exc_mass > 0:
                    exc_pts = int((exc_mass / comp_total_mass) * total_particles_per_comp)
                    if exc_pts > 0:
                        exc_coords = organ_coords['excretion']
                        px = np.random.uniform(exc_coords[0] - exc_coords[2]/2, exc_coords[0] + exc_coords[2]/2, exc_pts)
                        py = np.random.uniform(exc_coords[1] - exc_coords[3]/2, exc_coords[1] + exc_coords[3]/2, exc_pts)
                        particle_x.extend(px); particle_y.extend(py)

                if c_idx == 0 and git_mass > 0:
                    git_pts = int((git_mass / comp_total_mass) * total_particles_per_comp)
                    if git_pts > 0:
                        theta = np.random.uniform(0, 2*np.pi, git_pts)
                        r = np.sqrt(np.random.uniform(0, 1, git_pts))
                        px_git = 0.0 + r * 0.3 * np.cos(theta)
                        py_git = 5.15 + r * 1.65 * np.sin(theta)
                        particle_x.extend(px_git); particle_y.extend(py_git)

                scatter_handles[c_idx].set_offsets(np.c_[particle_x, particle_y])
                scatter_handles[c_idx].set_alpha(max(0.2, min(1.0, comp_total_mass / (dose_mg * 0.8))) * 0.5)

            ax_body.set_title(f'Simulation: {drug_name} (Dosing + Tail)\nTime: {sol_t[frame]:.1f}h | Total Mass: {total_body_mass:.2f} mg', 
                         color='black', fontsize=11, fontweight='bold', pad=12)
            return scatter_handles

        ani_body = animation.FuncAnimation(fig_body, update_body, frames=len(sol_t), interval=100, repeat=False)
        animations.append(('body', ani_body))

    if anim_type in ['pk', 'both']:
        fig_pk, ax_pk = plt.subplots(figsize=(10, 6))
        fig_pk.subplots_adjust(right=0.7) 
        
        full_plasma_conc = sol_y[venous_idx, :] / (params['V_blood'] * 0.67 * params['BP_Ratio'])
        max_y = c_max
        min_y = max(target_threshold / 5, 1e-4)

        lines_pk = []
        line_p, = ax_pk.plot([], [], label=f'Venous Plasma (Parent: {drug_name})', color='#0088cc', lw=2)
        lines_pk.append((line_p, np.maximum(full_plasma_conc, 1e-6)))

        met_colors = ['#ff6600', '#009966', '#cc00cc', '#d62728'] 
        for m_idx, m_params in enumerate(met_params_list):
            m_venous_idx = 1 + n_tissues * (m_idx + 1) + tissues.index('venous_blood')
            m_plasma_conc = sol_y[m_venous_idx, :] / (params['V_blood'] * 0.67 * m_params['BP_Ratio'])
            max_y = max(max_y, np.max(m_plasma_conc))
            
            m_name = m_params['Drug']
            c = met_colors[m_idx % len(met_colors)]
            line_m, = ax_pk.plot([], [], label=f'Venous Plasma ({m_name})', color=c, lw=2, linestyle='-')
            lines_pk.append((line_m, np.maximum(m_plasma_conc, 1e-6)))

        ax_pk.axhline(target_threshold, color='red', linestyle='--', label='1/100th Parent Cmax Threshold')
        ax_pk.set_yscale('log')
        ax_pk.set_ylim(bottom=min_y, top=max_y * 3)
        ax_pk.set_xlim(0, sol_t[-1] * 1.02)
        ax_pk.set_xlabel('Time (hours)')
        ax_pk.set_ylabel('Concentration (mg/L)')
        ax_pk.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        ax_pk.grid(True, alpha=0.3, which="both", ls="--")

        def update_pk(frame):
            for line_obj, data_y in lines_pk:
                line_obj.set_data(sol_t[:frame+1], data_y[:frame+1])
            ax_pk.set_title(f'PK Profile Animation: {drug_name} (Dosing + Tail)\nTime: {sol_t[frame]:.1f} hours')
            return [l[0] for l in lines_pk]

        ani_pk = animation.FuncAnimation(fig_pk, update_pk, frames=len(sol_t), interval=100, repeat=False)
        animations.append(('pk', ani_pk))

    global global_anis
    global_anis = [a[1] for a in animations]

    if save_gif:
        for name, ani in animations:
            gif_name = f"{drug_name}_{name}_simulation.gif"
            ani.save(gif_name, writer='pillow', fps=15)
            print(f"[+] Animation saved successfully as {gif_name}")
    else:
        if animations:
            plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Anatomical PBPK ADME Animator")
    parser.add_argument('--drug', type=str, default='Efavirenz', help='Parent drug name')
    parser.add_argument('--dose', type=float, default=100.0, help='Dose in mg')
    parser.add_argument('--duration_days', type=float, default=30.0, help='Dosing duration in days')
    parser.add_argument('--interval_hours', type=float, default=24.0, help='Dosing interval in hours')
    parser.add_argument('--anim_type', type=str, choices=['body', 'pk', 'both'], default='body', help='Which animation to run')
    parser.add_argument('--save_gif', action='store_true', help='Save animation as GIF')
    parser.add_argument('--drug_csv', type=str, default='benchmark_drugs.csv', help='Drug CSV database path')
    parser.add_argument('--met_csv', type=str, default='metabolites.csv', help='Metabolite CSV database path')
    args = parser.parse_args()
    
    run_animator(args.drug, args.dose, args.duration_days, args.interval_hours, args.anim_type, args.save_gif, args.drug_csv, args.met_csv)