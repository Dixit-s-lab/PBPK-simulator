PBPK-QSP Simulation Suite

An advanced Physiologically Based Pharmacokinetic (PBPK) and Quantitative Systems Pharmacology (QSP) modeling suite. 

This tool integrates a 12-state solid/dissolved ACAT gut model, Michaelis-Menten active transporter kinetics, Enterohepatic Recirculation (EHC), non-linear population scaling, and dynamic mid-simulation physiological events.
Features Thermodynamic Solubility: 

Auto-calculates exact solubility limits ($C_s$) from SMILES using RDKit's ESOL algorithm and applies Noyes-Whitney dissolution constraints.

Active Transport & EHC: Mechanistic tracking of biliary efflux and active renal secretion using saturable $V_{max}$ and $K_m$ kinetics.

Population & Disease Scaling: Dynamically scales organ volumes (e.g., obesity) and modifies clearances for Cirrhosis (Mild/Mod/Severe) and CKD (Stage 3/Stage 4).

QSP Event Scheduler: Triggers mid-simulation physiological changes (e.g., High-Fat Meals, DDI Inhibitor/Inducer start/stop, Smoking cessation).

Automated NCA & Validation: Calculates $C_{max}$, $T_{max}$, AUC, $t_{1/2}$, $F$, exact biological clearances, and generates $RMSE$/$R^2$ metrics against uploaded experimental data.

Split-Screen Dashboard: A fully interactive Tkinter GUI for parameter management and real-time console tracking.

🛠 InstallationBecause this suite relies heavily on RDKit for cheminformatics and ESOL solubility predictions, using a Conda environment is highly recommended to prevent C++ dependency conflicts.

1. Clone or create the directory:
2. Create a folder named poppbpk and place all 5 Python/CSV files inside it.
3. Create a fresh Conda environment:Bashconda create -n pbpk_env python=3.10

conda activate pbpk_env

4. Install RDKit:Bash
   conda install -c conda-forge rdkit
5. Install remaining Python dependencies:Bash
   pip install numpy pandas scipy matplotlib pubchempy openpyxl
   
📁 Repository Structure

gui1.py — The graphical user interface and primary entry point.

master_controllerACAT.py — The core ODE mathematical engine, NCA analyzer, and data exporter.

animatorACAT.py — The real-time Matplotlib Tkinter visualization engine for body flux and PK curves.

physiology_builder.py — Handles non-linear allometric population scaling and disease state modifiers.

distribution.py — Calculates tissue-to-plasma partition coefficients (Kp) via the Poulin & Theil mechanistic model.

patient_events.csv — The timeline scheduler for dynamic QSP events.


🚀 Usage & Sample CasesTo start the software, activate your environment and launch the dashboard:Bash
python gui1.py

Sample Case 1: Baseline PK Prediction (Auto-Fetch)Test how the software dynamically fetches drug properties from PubChem and runs a standard baseline simulation.

In the GUI Dosing Parameters, set Drug Name to Cetirizine.
Set Dose to 10.0, Duration to 3.0 days, and Interval to 24.0 hours.
Leave Patient Covariates at standard (70kg, None, Healthy).
Click Run Master Controller.

Result: The console will notify you that it fetched Cetirizine's SMILES, calculated its ESOL limit, generated a default metabolite, and saved it to benchmark_drugs.csv. It will output the full NCA analysis in the console and generate Cetirizine_Phase4_Plasma_Conc.csv and a PNG plot in your directory.

Sample Case 2: Population & Disease Scaling (PopPK)Test how liver disease and obesity alter drug accumulation and systemic clearance.
In the GUI, set Drug Name to Efavirenz. Set Dose to 600.0.
Under Patient Covariates, change Weight to 120.0 kg and Cirrhosis to Severe.
Click Run Master Controller.Result: You will see the total Hepatic Mass Metabolized drop significantly compared to a baseline patient, leading to a massive increase in $C_{max}$ and AUC. Tissue accumulation in the Adipose compartment will scale non-linearly due to the weight adjustment.

Sample Case 3: Dynamic QSP Event (Mid-Simulation DDI)Test the software's ability to trigger real-time physiological shifts. Open patient_events.csv in a text editor or Excel and add the following rows:

Time_hr,Event_Type,Value
48.0,DDI_Inhibitor,Started
96.0,DDI_Inhibitor,Stopped

In the GUI, set the Drug to Efavirenz, Duration to 7.0 days, and Dose to 100.0. Ensure the Event Timeline entry points to patient_events.csv. 
Click Run Animator. 

Select "Generate Both Animations Simultaneously" in the options.Result: Watch the animation dashboard. 
For the first two days, metabolism is normal. Exactly at Hour 48, an event will trigger in the console: Hepatic clearance crashes by 80%. The parent drug concentration will suddenly spike, and metabolite formation will stall. At Hour 96, the inhibitor stops, and the liver rapidly clears the backlog.

Sample Case 4: Experimental Data Validation ($R^2$ / $RMSE$)Test the mathematical correlation engine against in vivo clinical data.Create a simple validation file named clinical_data.csv in your folder. Format it with Time (hours) in the first column and Concentration (mg/L) in the second:

Time,Concentration
1.0,0.52
2.0,1.15
4.0,2.30
8.0,1.80
24.0,0.45

In the GUI Data Sources, use the Browse button next to Experimental File to select clinical_data.csv.
Click Run Master Controller.

Result: The system will linearly interpolate the simulated ODE curve to match your exact experimental time points. It will print the $R^2$ correlation and $RMSE$ to the console, and overlay the experimental scatter points as black dots directly onto the exported PNG PK curve.

⚖️ Licensing & UsageFor Academia & Students:This software is provided free of charge under a Custom Proprietary License strictly for non-commercial, educational, and academic thesis purposes.For Private Companies & Commercial Use:Commercial use, internal corporate R&D, software modification, and distribution are strictly prohibited without a separate commercial agreement. For commercial licensing and pricing, please contact the repository owner.
