# SEPIC vs ZETA Converter Analysis

This repo contains the MATLAB, PSIM, and Python files used to compare a 12V / 50W SEPIC and ZETA converter running at 340 kHz. 

The goal was to build analytical loss models in MATLAB, simulate them in PSIM, and use Python to see how well the math matches the waveforms.

## What's in here

- **`matlab/`** — The analytical models. These calculate component losses, efficiency, and junction temperatures over nine operating points. Also includes the Monte Carlo tolerance sweep and the thermal transient script.
- **`psim/`** — The PSIM schematics and run scripts. These run the steady-state simulations and the double-pulse test (DPT) to extract switching energies.
- **`python/`** — The validation scripts. These read the MATLAB CSVs and the PSIM exports, merge them, and spit out the error figures and summary tables.

## What you need before running

- **MATLAB R2020a+** (Needs `tiledlayout` and `exportgraphics`). You'll also need the Statistics toolbox for the Monte Carlo script (`lhsdesign`).
- **PSIM 2021+**
- **Python 3.9+** with `pandas`, `numpy`, and `matplotlib`.

## The Workflow (Run it in this order)

If you're running this from scratch, do it in this order or the Python scripts will crash looking for files that don't exist yet.

### 1. Generate the DPT Lookup Table
The MATLAB models rely on `dpt_lookup.csv`, which is the measured switching energies from the double-pulse test. 
1. Open `psim/DPT.psimsch` and run its sweep script. 
2. It outputs a CSV. Rename it to `dpt_lookup.csv` and put it in the working directory when you run MATLAB. If you skip this, the MATLAB scripts will fall back to a dummy 28V/5A point and your switching losses will be wrong.

### 2. Run the PSIM Steady-State Sweeps
1. Open the `Clear();.txt` script in `psim/SEPIC NEW/` and `psim/ZETA NEW/` from the PSIM script window.
2. Check the `SchFile` path at the top. If the file on your disk is called `SEPIC.psimsch` but the script says `FINAL SEPIC .psimsch` (watch the space!), it will fail.
3. Run it. This creates the `.txt` and `.smv` files for all 9 cases.

*Note: The TXT exports are at 2ns timestep so the Python script can integrate the switching loss correctly. The SMV files are at 10ns just so you can look at them in the PSIM viewer without your PC freezing.*

### 3. Run the MATLAB Models
1. Open `matlab/SEPIC_ANALYTICAL_MODEL_v9_16.m` or the ZETA equivalent.
2. Make sure `dpt_lookup.csv` is in your MATLAB path.
3. Run it. It will generate the summary tables and the IEEE/report figures in `report_figures/`.

### 4. Run the Python Validation
Once both PSIM and MATLAB have generated their CSVs, run the Python scripts in `python/`.
1. Edit the file paths at the top of the scripts so `MATLAB_DIR` and `PSIM_OUTPUT_DIR` point to your actual folders.
2. Run `SEPIC_CORRECTED_PSIM_ANALYSIS.py` to process the raw PSIM text files.
3. Run the validation script (e.g., `PSIM VS MATLAB VALIDATION SCRIPT.py`). It will generate the V1-V4 comparison figures.

## Watch out for these things

- **The PSIM switching window:** The Python script looks at a 200ns/500ns window around the switching edges to calculate loss. ZETA's PWM-to-transition delay is longer than SEPIC's, so if you see switching loss coming out near zero, the window missed the crossover event. You can widen `SWITCHING_POST_NS` in the Python script if needed.
- **DPT CSV Columns:** The MATLAB loader expects specific column names (`I_turnoff_A`, `I_turnon_A`, etc.). If your DPT export doesn't have them, it falls back to a scalar reference current. 
- **Font errors:** If you run the MATLAB scripts on a Linux machine without Times New Roman, the `ieee_cfg()` function will fall back to whatever serif font it can find. The figures will still generate, they just won't be IEEE-strict.

## License
MIT. Use it if it's helpful.
