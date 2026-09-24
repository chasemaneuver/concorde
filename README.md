# Analysis and Resizing of the Concorde Rolls-Royce/Snecma Olympus 593MK610

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22068665.svg)](https://doi.org/10.5281/zenodo.22068665)

This repository contains the MATLAB and Python codes developed for the academic project *Analysis and Resizing of the Concorde Rolls-Royce/Snecma Olympus 593MK610*, carried out as part of the final coursework requirements for Aerospace Engineering at Politecnico di Milano.

The project investigates the Rolls-Royce/Snecma Olympus 593MK610 turbojet installed on the Aerospatiale-BAC Concorde. The original engine cycle is analyzed under four representative operating conditions: take-off, noise-abatement climb at sea level, maximum climb at an altitude of 12,000 m, and supersonic cruise. The study covers the external-compression intake and its shock-wave system, compressor velocity triangles, engine-cycle thermodynamics, coaxial exhaust flow, nozzle operation, thrust, specific fuel consumption, specific impulse, and component and overall efficiencies.

A preliminary redesign is also examined. The proposed configuration converts the Olympus into a low-bypass, associated-flow turbofan while retaining a tertiary stream for excess-flow management, cooling, mixing, and noise control. A turboramjet derivative and chevron-equipped nozzle concepts are additionally considered for comparison.

The variable names and symbols used in the scripts follow, wherever possible, the notation adopted in the report. The codes are provided to document and support the calculations, assumptions, and numerical results presented in the associated publication.

## Analysis scope

The computational work includes:

- thermodynamic cycle calculations for the original Olympus 593MK610;
- performance evaluation at take-off, noise-abatement climb, maximum climb, and supersonic cruise;
- oblique- and normal-shock calculations for the external-compression intake;
- total-pressure and total-temperature calculations across the main engine components;
- compressor and fan mean-line velocity-triangle analysis;
- preliminary rotor-blade sizing and material-stress calculations;
- convergent and convergent-divergent nozzle calculations;
- thrust, specific impulse, TSFC, and efficiency evaluation;
- analysis of low-bypass turbofan and turboramjet derivatives;
- parametric investigation of bypass-ratio and turbine-entry-temperature effects.

GasTurb13, SolidWorks, and ANSYS were also used in the wider project for engine-cycle simulation, geometric design, and CFD analysis. Their proprietary model files are not included in this repository.

## Numerical approach and limitations

The scripts implement preliminary one-dimensional or quasi-one-dimensional engine-cycle and component models. They use steady operating conditions, ideal-gas relations with assigned thermodynamic properties, prescribed component efficiencies and pressure ratios, and the simplifying assumptions documented in the report and source-code comments.

The compressor and fan calculations rely on mean-line velocity triangles and repeated-stage or repeated-row assumptions where indicated. The intake calculations use classical compressible-flow relations for oblique and normal shocks, while the nozzle calculations use critical-flow and adapted-flow conditions.

The results are intended for academic comparison and preliminary design. They do not constitute a certified engine-performance model, detailed off-design analysis, structural substantiation, or complete propulsion-system redesign.

## Selected results

### Original Olympus 593MK610

| Operating condition | Flight Mach number | Thrust (kN) | Specific impulse (s) | TSFC (x 10^-5 kg/(N s)) | Overall efficiency |
|---|---:|---:|---:|---:|---:|
| Supersonic cruise | 2.0000 | 44.6670 | 2819 | 36.1600 | 0.3783 |
| Take-off | 0.3320 | 150.3800 | 3101 | 32.8000 | 0.0797 |
| Noise-abatement climb at sea level | 0.3775 | 173.0700 | 3093 | 32.9000 | 0.0904 |
| Maximum climb at 12,000 m | 1.2000 | 91.2810 | 2549 | 39.9000 | 0.2053 |

### Supersonic-cruise configuration comparison

| Configuration | Flight Mach number | Thrust (kN) | Specific impulse (s) | TSFC (x 10^-5 kg/(N s)) | Thermal efficiency | Propulsive efficiency | Overall efficiency |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original Olympus | 2 | 44.6670 | 2819 | 36.1600 | 0.5157 | 0.7335 | 0.3783 |
| Proposed turbofan | 2 | 70.5480 | 3831 | 26.6000 | 0.7709 | 0.6669 | 0.5141 |
| Turboramjet derivative | 3 | 45.9060 | 1957 | 52.1000 | 0.5053 | 0.7796 | 0.3940 |

Within the assumptions of the preliminary model, the proposed turbofan configuration increases thrust and overall efficiency while reducing TSFC relative to the original engine at the analyzed cruise condition. The report also presents a preliminary CFD study of chevron-equipped nozzles, indicating a noise reduction of approximately 10.8 dB at the cost of a small thrust penalty.

## Repository structure

```text
.
|-- matlab/
|   |-- Cruise.m
|   |-- Takeoff.m
|   |-- Noise_abatement_climb.m
|   |-- Max_climb.m
|   |-- Derive_turbofan.m
|   |-- Derived_turboramjet.m
|   `-- function/
|       |-- Inlet_and_shockwaves.m
|       |-- Nozzle.m
|       `-- Nozzle_inlet_critical_flow.m
|-- python/
|   |-- supersonic_cruise.py
|   |-- compressor_analysis.py
|   |-- fan_design.py
|   |-- rotor_blade_design.py
|   |-- parameter_misc.py
|   |-- Test.py
|   `-- function/
|       |-- intro_parameters.py
|       |-- parameters_engine_components.py
|       `-- shock_waves.py
|-- LICENSE
`-- README.md
```

### MATLAB files

- `Cruise.m`, `Takeoff.m`, `Noise_abatement_climb.m`, and `Max_climb.m` evaluate the original engine at the four operating conditions examined in the report.
- `Derive_turbofan.m` evaluates the proposed associated-flow turbofan derivative.
- `Derived_turboramjet.m` evaluates the preliminary turboramjet derivative.
- The `matlab/function` directory contains the intake, shock-wave, nozzle, and critical-flow helper functions.

### Python files

- `supersonic_cruise.py` reproduces the supersonic-cruise engine-cycle analysis.
- `compressor_analysis.py` and `fan_design.py` implement the mean-line turbomachinery analyses.
- `rotor_blade_design.py` supports preliminary rotor-blade sizing and material calculations.
- `parameter_misc.py` contains thrust, TSFC, specific-impulse, and efficiency utilities.
- `Test.py` contains exploratory calculations for the proposed turbofan configuration.
- The `python/function` directory contains atmospheric, engine-component, and shock-wave functions.

## Software requirements

- MATLAB;
- Python 3;
- NumPy;
- Matplotlib.

Install the required Python packages with:

```bash
python -m pip install numpy matplotlib
```

GasTurb13, SolidWorks, and ANSYS are only required to reproduce the corresponding simulations and design work described in the report; they are not required for the MATLAB and Python scripts stored here.

## Important setup notes

The repository preserves several filenames and import names from the original academic project. Before running the scripts, make the following names consistent.

### MATLAB helper functions

MATLAB requires the primary function name to match its filename. Either rename the files as follows or update the corresponding function declarations and calls:

| Current filename | Function name used by the scripts | Suggested filename |
|---|---|---|
| `Inlet_and_shockwaves.m` | `inlet10` | `inlet10.m` |
| `Nozzle.m` | `nozzle10` | `nozzle10.m` |
| `Nozzle_inlet_critical_flow.m` | `nozzleinlet` | `nozzleinlet.m` |

Add the MATLAB folders to the path before execution:

```matlab
addpath(genpath('matlab'));
```

### Python helper modules

The Python scripts use compact legacy module names that do not exactly match several current filenames. Either update the import statements to use the filenames shown in the repository or rename the modules consistently. The main correspondences are:

| Current filename | Legacy import name |
|---|---|
| `intro_parameters.py` | `introparameters` |
| `parameters_engine_components.py` | `parametersenginecomponents` |
| `shock_waves.py` | `shockwaves` |
| `parameter_misc.py` | `parametersmisc` |
| `compressor_analysis.py` | `compressoranalysis` |

The `python/function` directory must also be included in the Python module search path.

## Usage

Clone or download the repository and preserve its folder structure.

For the MATLAB analyses, open the repository root in MATLAB, add the `matlab` directory and its subdirectories to the path, apply the helper-file naming adjustments described above, and run the script associated with the desired operating condition or engine configuration.

For the Python analyses, make the module names and imports consistent, ensure that both `python` and `python/function` are accessible, and execute the desired script. For example:

```bash
python python/supersonic_cruise.py
```

Before changing an input case, verify the units and assumptions used by the selected script. The original files contain a mixture of SI units and engineering-unit scalings consistent with the tables in the report.

## Associated publication

The complete report is available on Zenodo:

**DOI:** [10.5281/zenodo.22068665](https://doi.org/10.5281/zenodo.22068665)

## Citation

If you use the codes or results contained in this repository, please cite the associated report:

> M. Vrapi, A. R. Shah, P. Roncoroni, A. Perego, and B. Tassinari, "Analysis and Resizing of the Concorde Rolls-Royce/Snecma Olympus 593MK610," Zenodo, 2021, doi: [10.5281/zenodo.22068665](https://doi.org/10.5281/zenodo.22068665).

## License

The source code in this repository is distributed under the MIT License. See the `LICENSE` file for the complete terms. The associated report is distributed separately under the Creative Commons Attribution 4.0 International License (CC BY 4.0), as specified in its Zenodo record.

## Authors

- Michelle Vrapi ([ORCID](https://orcid.org/0009-0007-3304-8042))
- Areeb Raza Shah
- Paolo Roncoroni
- Andrea Perego
- Beniamino Tassinari

Politecnico di Milano, Department of Aerospace Science and Technology (DAER).

## Project context

This repository is associated with an academic aerospace engineering project and is not an official publication of Politecnico di Milano. The codes are shared for documentation, reproducibility, and portfolio purposes.
