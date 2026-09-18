# gui1.py
import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog
import subprocess
import sys
import threading

class StdRedirector:
    def __init__(self, widget):
        self.widget = widget

    def write(self, str_data):
        self.widget.insert(tk.END, str_data)
        self.widget.see(tk.END)

    def flush(self):
        pass

class PBPKControlCenter:
    def __init__(self, root):
        self.root = root
        self.root.title("PBPK QSP Control Center (Phase 4)")
        # Expanded geometry to fit the side-by-side layout comfortably
        self.root.geometry("1300x800")
        self.root.configure(bg="#f4f6f9")

        title_label = tk.Label(root, text="PBPK Population & Event Simulation Suite", font=("Arial", 18, "bold"), bg="#f4f6f9", fg="#333")
        title_label.pack(pady=10)

        # --- Create a PanedWindow for Left/Right Split ---
        self.main_pane = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.main_pane.pack(fill="both", expand=True, padx=15, pady=10)

        # Left Panel (Inputs & Options)
        self.left_panel = tk.Frame(self.main_pane, bg="#f4f6f9")
        # Right Panel (Console & Execution)
        self.right_panel = tk.Frame(self.main_pane, bg="#f4f6f9")

        self.main_pane.add(self.left_panel, weight=1)
        self.main_pane.add(self.right_panel, weight=1)

        # =============================================================
        # LEFT PANEL: INPUTS & SETTINGS
        # =============================================================
        
        # --- Data Sources Frame ---
        data_frame = ttk.LabelFrame(self.left_panel, text="Data Sources (CSV)", padding=10)
        data_frame.pack(fill="x", pady=(0, 10))
        
        ttk.Label(data_frame, text="Drug DB:", font=("Arial", 10)).grid(row=0, column=0, sticky="w", pady=5)
        self.drug_csv_entry = ttk.Entry(data_frame, width=35, font=("Arial", 10))
        self.drug_csv_entry.insert(0, "benchmark_drugs.csv")
        self.drug_csv_entry.grid(row=0, column=1, sticky="w", pady=5, padx=5)
        ttk.Button(data_frame, text="Browse", command=self.browse_drug_csv).grid(row=0, column=2, padx=5)

        ttk.Label(data_frame, text="Metabolite DB:", font=("Arial", 10)).grid(row=1, column=0, sticky="w", pady=5)
        self.met_csv_entry = ttk.Entry(data_frame, width=35, font=("Arial", 10))
        self.met_csv_entry.insert(0, "metabolites.csv")
        self.met_csv_entry.grid(row=1, column=1, sticky="w", pady=5, padx=5)
        ttk.Button(data_frame, text="Browse", command=self.browse_met_csv).grid(row=1, column=2, padx=5)

        ttk.Label(data_frame, text="Event Timeline:", font=("Arial", 10)).grid(row=2, column=0, sticky="w", pady=5)
        self.events_csv_entry = ttk.Entry(data_frame, width=35, font=("Arial", 10))
        self.events_csv_entry.insert(0, "patient_events.csv")
        self.events_csv_entry.grid(row=2, column=1, sticky="w", pady=5, padx=5)
        ttk.Button(data_frame, text="Browse", command=self.browse_events_csv).grid(row=2, column=2, padx=5)

        # --- Dosing Frame ---
        input_frame = ttk.LabelFrame(self.left_panel, text="Dosing Parameters", padding=10)
        input_frame.pack(fill="x", pady=10)

        ttk.Label(input_frame, text="Drug Name:", font=("Arial", 10)).grid(row=0, column=0, sticky="w", pady=5)
        self.drug_entry = ttk.Entry(input_frame, width=15, font=("Arial", 10))
        self.drug_entry.insert(0, "Efavirenz")
        self.drug_entry.grid(row=0, column=1, sticky="w", pady=5, padx=5)

        ttk.Label(input_frame, text="Dose (mg):", font=("Arial", 10)).grid(row=0, column=2, sticky="w", pady=5)
        self.dose_entry = ttk.Entry(input_frame, width=15, font=("Arial", 10))
        self.dose_entry.insert(0, "100.0")
        self.dose_entry.grid(row=0, column=3, sticky="w", pady=5, padx=5)

        ttk.Label(input_frame, text="Duration (Days):", font=("Arial", 10)).grid(row=1, column=0, sticky="w", pady=5)
        self.duration_entry = ttk.Entry(input_frame, width=15, font=("Arial", 10))
        self.duration_entry.insert(0, "30.0")
        self.duration_entry.grid(row=1, column=1, sticky="w", pady=5, padx=5)

        ttk.Label(input_frame, text="Interval (Hours):", font=("Arial", 10)).grid(row=1, column=2, sticky="w", pady=5)
        self.interval_entry = ttk.Entry(input_frame, width=15, font=("Arial", 10))
        self.interval_entry.insert(0, "24.0")
        self.interval_entry.grid(row=1, column=3, sticky="w", pady=5, padx=5)

        # --- Patient Profile Frame ---
        patient_frame = ttk.LabelFrame(self.left_panel, text="Patient Covariates (Disease Scaling)", padding=10)
        patient_frame.pack(fill="x", pady=10)

        ttk.Label(patient_frame, text="Weight (kg):", font=("Arial", 10)).grid(row=0, column=0, sticky="w", pady=5)
        self.weight_entry = ttk.Entry(patient_frame, width=15, font=("Arial", 10))
        self.weight_entry.insert(0, "70.0")
        self.weight_entry.grid(row=0, column=1, sticky="w", pady=5, padx=5)

        ttk.Label(patient_frame, text="Cirrhosis:", font=("Arial", 10)).grid(row=0, column=2, sticky="w", pady=5)
        self.cirrhosis_var = tk.StringVar(value="None")
        cirrhosis_dropdown = ttk.Combobox(patient_frame, textvariable=self.cirrhosis_var, values=["None", "Mild", "Moderate", "Severe"], state="readonly", width=12)
        cirrhosis_dropdown.grid(row=0, column=3, sticky="w", pady=5, padx=5)

        ttk.Label(patient_frame, text="Renal Status:", font=("Arial", 10)).grid(row=1, column=2, sticky="w", pady=5)
        self.renal_var = tk.StringVar(value="Healthy")
        renal_dropdown = ttk.Combobox(patient_frame, textvariable=self.renal_var, values=["Healthy", "CKD_Stage_3", "CKD_Stage_4"], state="readonly", width=12)
        renal_dropdown.grid(row=1, column=3, sticky="w", pady=5, padx=5)

        # --- Options Frame ---
        options_frame = ttk.LabelFrame(self.left_panel, text="Execution Options", padding=10)
        options_frame.pack(fill="x", pady=10)

        self.save_png_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Save PNG Summary Plot (Master Controller)", variable=self.save_png_var).pack(anchor="w", pady=(2, 10))

        ttk.Label(options_frame, text="Select Animation Style (For Animator):", font=("Arial", 10, "bold")).pack(anchor="w")
        self.anim_type_var = tk.StringVar(value="both")
        ttk.Radiobutton(options_frame, text="Body Organ Level Only", variable=self.anim_type_var, value="body").pack(anchor="w")
        ttk.Radiobutton(options_frame, text="PK Concentration Time Curve Only", variable=self.anim_type_var, value="pk").pack(anchor="w")
        ttk.Radiobutton(options_frame, text="Generate Both Animations Simultaneously", variable=self.anim_type_var, value="both").pack(anchor="w")

        self.save_gif_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Save Animation(s) as GIF instead of displaying", variable=self.save_gif_var).pack(anchor="w", pady=(10, 2))

        # =============================================================
        # RIGHT PANEL: EXECUTION BUTTONS & CONSOLE
        # =============================================================
        
        # --- Buttons Frame ---
        btn_frame = tk.Frame(self.right_panel, bg="#f4f6f9")
        btn_frame.pack(fill="x", pady=(0, 10))

        # Making buttons prominent
        run_mc_btn = ttk.Button(btn_frame, text="Run Master Controller", command=self.run_master_controller)
        run_mc_btn.pack(side="left", padx=2, expand=True, fill="both", ipady=8)

        run_anim_btn = ttk.Button(btn_frame, text="Run Animator", command=self.run_animator)
        run_anim_btn.pack(side="left", padx=2, expand=True, fill="both", ipady=8)

        clear_btn = ttk.Button(btn_frame, text="Clear Console", command=self.clear_console)
        clear_btn.pack(side="left", padx=2, expand=True, fill="both", ipady=8)

        # --- Execution Console ---
        console_frame = ttk.LabelFrame(self.right_panel, text="Execution Console & Log", padding=5)
        console_frame.pack(fill="both", expand=True)

        self.console = scrolledtext.ScrolledText(console_frame, wrap=tk.WORD, bg="#1e1e1e", fg="#00ff00", font=("Courier", 11))
        self.console.pack(fill="both", expand=True)

        self.redirector = StdRedirector(self.console)

    def browse_drug_csv(self):
        filename = filedialog.askopenfilename(title="Select Drug Database CSV", filetypes=[("CSV files", "*.csv")])
        if filename:
            self.drug_csv_entry.delete(0, tk.END)
            self.drug_csv_entry.insert(0, filename)

    def browse_met_csv(self):
        filename = filedialog.askopenfilename(title="Select Metabolite Database CSV", filetypes=[("CSV files", "*.csv")])
        if filename:
            self.met_csv_entry.delete(0, tk.END)
            self.met_csv_entry.insert(0, filename)

    def browse_events_csv(self):
        filename = filedialog.askopenfilename(title="Select Events Timeline CSV", filetypes=[("CSV files", "*.csv")])
        if filename:
            self.events_csv_entry.delete(0, tk.END)
            self.events_csv_entry.insert(0, filename)

    def clear_console(self):
        self.console.delete('1.0', tk.END)

    def run_master_controller(self):
        drug = self.drug_entry.get().strip()
        dose = self.dose_entry.get().strip()
        duration = self.duration_entry.get().strip()
        interval = self.interval_entry.get().strip()
        weight = self.weight_entry.get().strip()
        cirrhosis = self.cirrhosis_var.get()
        renal = self.renal_var.get()
        
        save_png = self.save_png_var.get()
        drug_csv = self.drug_csv_entry.get().strip()
        met_csv = self.met_csv_entry.get().strip()
        events_csv = self.events_csv_entry.get().strip()

        def task():
            cmd = [sys.executable, "master_controllerACAT.py", 
                   "--drug", drug, "--dose", dose, 
                   "--duration_days", duration, "--interval_hours", interval, 
                   "--weight", weight, "--cirrhosis", cirrhosis, "--renal", renal,
                   "--drug_csv", drug_csv, "--met_csv", met_csv, "--events_csv", events_csv]
            if save_png:
                cmd.append("--save_png")
            
            self.console.insert(tk.END, f"\n>>> Executing: {' '.join(cmd)}\n")
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                self.redirector.write(line)
            process.wait()
            self.console.insert(tk.END, "\n[+] Master Controller execution finished.\n")

        threading.Thread(target=task, daemon=True).start()

    def run_animator(self):
        drug = self.drug_entry.get().strip()
        dose = self.dose_entry.get().strip()
        duration = self.duration_entry.get().strip()
        interval = self.interval_entry.get().strip()
        weight = self.weight_entry.get().strip()
        cirrhosis = self.cirrhosis_var.get()
        renal = self.renal_var.get()
        
        anim_type = self.anim_type_var.get()
        save_gif = self.save_gif_var.get()
        
        drug_csv = self.drug_csv_entry.get().strip()
        met_csv = self.met_csv_entry.get().strip()
        events_csv = self.events_csv_entry.get().strip()

        def task():
            cmd = [sys.executable, "animatorACAT.py", 
                   "--drug", drug, "--dose", dose, 
                   "--duration_days", duration, "--interval_hours", interval, 
                   "--weight", weight, "--cirrhosis", cirrhosis, "--renal", renal,
                   "--anim_type", anim_type, 
                   "--drug_csv", drug_csv, "--met_csv", met_csv, "--events_csv", events_csv]
            if save_gif:
                cmd.append("--save_gif")

            self.console.insert(tk.END, f"\n>>> Executing Animator ({anim_type} mode): {' '.join(cmd)}\n")
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                self.redirector.write(line)
            process.wait()
            self.console.insert(tk.END, "\n[+] Animator execution finished.\n")

        threading.Thread(target=task, daemon=True).start()

if __name__ == "__main__":
    root = tk.Tk()
    app = PBPKControlCenter(root)
    root.mainloop()