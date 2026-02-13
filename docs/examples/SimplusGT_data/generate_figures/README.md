# Figure Generation Scripts

This folder contains scripts to reproduce the analysis figures for the Simplus Grid Tool disturbance simulation study.

## Prerequisites

- MATLAB R2020a or later
- Simulink
- Git

## Setup Instructions

### 1. Download SimplusGT

Clone the SimplusGT repository at the specific commit used for this analysis:

```bash
git clone https://github.com/ANL-CEEESA/SimplusGT.git
cd SimplusGT
git checkout 30cf90c424ef116b69a43cafeefe4b26ff936886
```

### 2. Install SimplusGT

Run the installation script **once** in MATLAB:

```matlab
InstallSimplusGT
```

This will:
- Add SimplusGT folders to your MATLAB path
- Convert the Simulink library to your MATLAB version
- Set up the toolbox for use

### 3. Copy Figure Generation Scripts

Copy the `generate_figures` folder to the SimplusGT project root:

```
SimplusGT/
├── UserMain.m
├── InstallSimplusGT.m
├── Examples/
├── Library/
├── generate_figures/          # ← This folder
│   ├── README.md
│   ├── RunDisturbanceSim.m
│   └── ExportFigures.m
└── ...
```

## Generating Figures

Navigate to the SimplusGT root directory in MATLAB, then run the following scripts **in order**:

- `UserMain`
- `RunDisturbanceSim`
- `ExportFigures`
