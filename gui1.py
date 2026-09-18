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
        self.root.title("PBPK-VAD Integrated Control Center")
        self.root.geometry("760x780")
        self.root.configure(bg="#f4f6f9")

        title_label = tk.Label(root, text="PBPK ADME Simulation & Visualization Suite", font=("Arial", 16, "bold"), bg="#f4f6f9", fg="#333")
        title_label.pack(pady=10)

        # -------------------------------------------------------------
        # Data Sources Frame (NEW)
        # -------------------------------------------------------------
        data_frame = ttk.LabelFrame(root, text="Data Sources (CSV)", padding=15)
        data_frame.pack(fill="x", padx=20, pady=5)
        
        ttk.Label(data_frame, text="Drug DB:", font=("Arial", 10)).grid(row=0, column=0, sticky="w", pady=5)
        self.drug_csv_entry = ttk.Entry(data_frame, width=45, font=("Arial", 10))
        self.drug_csv_entry.insert(0, "benchmark_drugs.csv")
        self.drug_csv_entry.grid(row=0, column=1, sticky="w", pady=5, padx=5)
        ttk.Button(data_frame, text="Browse", command=self.browse_drug_csv).grid(row=0, column=2, padx=5)

        ttk.Label(data_frame, text="Metabolite DB:", font=("Arial", 10)).grid(row=1, column=0, sticky="w", pady=5)
        self.met_csv_entry = ttk.Entry(data_frame, width=45, font=("Arial", 10))
        self.met_csv_entry.insert(0, "metabolites.csv")
        self.met_csv_entry.grid(row=1, column=1, sticky="w", pady=5, padx=5)
        ttk.Button(data_frame, text="Browse", command=self.browse_met_csv).grid(row=1, column=2, padx=5)

        # -------------------------------------------------------------
        # Input Frame
        # -------------------------------------------------------------
        input_frame = ttk.LabelFrame(root, text="Simulation Parameters", padding=15)
        input_frame.pack(fill="x", padx=20, pady=5)

        ttk.Label(input_frame, text="Drug Name:", font=("Arial", 11)).grid(row=0, column=0, sticky="w", pady=5)
        self.drug_entry = ttk.Entry(input_frame, width=20, font=("Arial", 11))
        self.drug_entry.insert(0, "Efavirenz")
        self.drug_entry.grid(row=0, column=1, sticky="w", pady=5, padx=10)

        ttk.Label(input_frame, text="Dose (mg):", font=("Arial", 11)).grid(row=0, column=2, sticky="w", pady=5)
        self.dose_entry = ttk.Entry(input_frame, width=10, font=("Arial", 11))
        self.dose_entry.insert(0, "100.0")
        self.dose_entry.grid(row=0, column=3, sticky="w", pady=5, padx=10)

        ttk.Label(input_frame, text="Duration (Days):", font=("Arial", 11)).grid(row=1, column=0, sticky="w", pady=5)
        self.duration_entry = ttk.Entry(input_frame, width=20, font=("Arial", 11))
        self.duration_entry.insert(0, "30.0")
        self.duration_entry.grid(row=1, column=1, sticky="w", pady=5, padx=10)

        ttk.Label(input_frame, text="Interval (Hours):", font=("Arial", 11)).grid(row=1, column=2, sticky="w", pady=5)
        self.interval_entry = ttk.Entry(input_frame, width=10, font=("Arial", 11))
        self.interval_entry.insert(0, "24.0")
        self.interval_entry.grid(row=1, column=3, sticky="w", pady=5, padx=10)

        # -------------------------------------------------------------
        # Options Frame
        # -------------------------------------------------------------
        options_frame = ttk.LabelFrame(root, text="Execution Options", padding=15)
        options_frame.pack(fill="x", padx=20, pady=5)

        self.save_png_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Save PNG Summary Plot (Master Controller)", variable=self.save_png_var).pack(anchor="w", pady=(2, 10))

        ttk.Label(options_frame, text="Select Animation Style (For Animator):", font=("Arial", 10, "bold")).pack(anchor="w")
        self.anim_type_var = tk.StringVar(value="both")
        ttk.Radiobutton(options_frame, text="Body Organ Level Only", variable=self.anim_type_var, value="body").pack(anchor="w")
        ttk.Radiobutton(options_frame, text="PK Concentration Time Curve Only", variable=self.anim_type_var, value="pk").pack(anchor="w")
        ttk.Radiobutton(options_frame, text="Generate Both Animations Simultaneously", variable=self.anim_type_var, value="both").pack(anchor="w")

        self.save_gif_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Save Animation(s) as GIF instead of displaying", variable=self.save_gif_var).pack(anchor="w", pady=(10, 2))

        # -------------------------------------------------------------
        # Buttons & Output Frame
        # -------------------------------------------------------------
        btn_frame = tk.Frame(root, bg="#f4f6f9")
        btn_frame.pack(fill="x", padx=20, pady=10)

        run_mc_btn = ttk.Button(btn_frame, text="Run Master Controller", command=self.run_master_controller)
        run_mc_btn.pack(side="left", padx=5, expand=True, fill="x")

        run_anim_btn = ttk.Button(btn_frame, text="Run Animator", command=self.run_animator)
        run_anim_btn.pack(side="left", padx=5, expand=True, fill="x")

        clear_btn = ttk.Button(btn_frame, text="Clear Console", command=self.clear_console)
        clear_btn.pack(side="left", padx=5, expand=True, fill="x")

        console_frame = ttk.LabelFrame(root, text="Execution Console & Log", padding=10)
        console_frame.pack(fill="both", expand=True, padx=20, pady=5)

        self.console = scrolledtext.ScrolledText(console_frame, wrap=tk.WORD, bg="#1e1e1e", fg="#00ff00", font=("Courier", 10))
        self.console.pack(fill="both", expand=True)

        self.redirector = StdRedirector(self.console)

    def browse_drug_csv(self):
        filename = filedialog.asksaveasfilename(title="Select or Create Drug Database CSV", defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if filename:
            self.drug_csv_entry.delete(0, tk.END)
            self.drug_csv_entry.insert(0, filename)

    def browse_met_csv(self):
        filename = filedialog.asksaveasfilename(title="Select or Create Metabolite Database CSV", defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if filename:
            self.met_csv_entry.delete(0, tk.END)
            self.met_csv_entry.insert(0, filename)

    def clear_console(self):
        self.console.delete('1.0', tk.END)

    def run_master_controller(self):
        drug = self.drug_entry.get().strip()
        dose = self.dose_entry.get().strip()
        duration = self.duration_entry.get().strip()
        interval = self.interval_entry.get().strip()
        save_png = self.save_png_var.get()
        drug_csv = self.drug_csv_entry.get().strip()
        met_csv = self.met_csv_entry.get().strip()

        def task():
            cmd = [sys.executable, "master_controller1.py", "--drug", drug, "--dose", dose, 
                   "--duration_days", duration, "--interval_hours", interval, 
                   "--drug_csv", drug_csv, "--met_csv", met_csv]
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
        anim_type = self.anim_type_var.get()
        save_gif = self.save_gif_var.get()
        drug_csv = self.drug_csv_entry.get().strip()
        met_csv = self.met_csv_entry.get().strip()

        def task():
            cmd = [sys.executable, "animator1.py", "--drug", drug, "--dose", dose, 
                   "--duration_days", duration, "--interval_hours", interval, 
                   "--anim_type", anim_type, "--drug_csv", drug_csv, "--met_csv", met_csv]
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