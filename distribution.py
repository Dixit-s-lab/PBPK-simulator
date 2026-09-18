# distribution.py

def calculate_tissue_kp(drug_params):
    """
    Calculates tissue-to-plasma partition coefficients (Kp) 
    using a simplified tissue-composition-based Poulin & Theil model.
    """
    logp = drug_params.get('LogP', 1.0)
    fup = drug_params.get('fup', 1.0)
    
    # Convert LogP to raw partition coefficient
    p_ow = 10 ** logp
    
    # Tissue Compositions: [Fraction Water, Fraction Neutral Lipids, Fraction Phospholipids]
    tissues = {
        'Adipose': [0.15, 0.79, 0.005],
        'Bone':    [0.45, 0.05, 0.002],
        'Brain':   [0.75, 0.05, 0.050],
        'Heart':   [0.75, 0.015, 0.015],
        'Kidney':  [0.78, 0.02, 0.016],
        'Liver':   [0.72, 0.03, 0.025],
        'Lungs':   [0.78, 0.015, 0.015],
        'Muscle':  [0.76, 0.01, 0.007],
        'Skin':    [0.70, 0.03, 0.010]
    }
    
    kp_dict = {}
    for tissue, (v_wt, v_nlt, v_npt) in tissues.items():
        # Poulin & Theil baseline equation
        kp = fup * (v_wt + (v_nlt * p_ow) + (v_npt * (0.3 * p_ow + 0.7)))
        
        # Apply a scalar if present in the database, ensure it never hits absolute zero
        scalar = drug_params.get('Kp_Scalar', 1.0)
        kp_dict[tissue] = max(kp * scalar, 0.01)
        
    # Blood Kp is 1.0 by definition (ratio of plasma to plasma)
    kp_dict['Venous_blood'] = 1.0
    kp_dict['Arterial_blood'] = 1.0
    
    return kp_dict