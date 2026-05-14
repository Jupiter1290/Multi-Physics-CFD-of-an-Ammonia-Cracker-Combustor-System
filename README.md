# CFD Investigation of a Thermally Coupled Ammonia Decomposer–Combustor System

## Overview

This repository contains the numerical framework, MATLAB analysis, ANSYS Fluent User-Defined Functions (UDFs), and presentation material developed for a CFD investigation of a thermally coupled ammonia decomposer–combustor system.

The project combines:
- porous media modeling,
- finite-rate chemistry,
- reacting-flow CFD,
- conjugate heat transfer,
- and custom numerical transport modeling using UDFs.

---

# Table of Contents

1. Project Objectives  
2. Physical System  
3. Numerical Methodology  
4. Repository Structure  
5. MATLAB Scripts  
6. UDF-Based Coupling Strategy  
7. Zone-Specific Diffusivity Modeling  
8. Key Results  
9. Software Used  

---

# Project Objectives

- Reproduce and validate the ammonia decomposition model of Chein et al.
- Perform standalone CFD simulations of the decomposer and hydrogen–air combustor to evaluate the feasibility of thermal integration.
- Redesign and radially scale the decomposer to accommodate the higher mass flow rates required by the combustor.
- Develop an annular decomposer configuration for flue-gas recirculation and thermal coupling.
- Perform standalone thermally coupled simulations by computing recirculated flue-gas inlet conditions using MATLAB and area-averaged combustor outlet properties.
- Implement fully coupled CFD simulations using UDFs for:
  - dynamic flue-gas recirculation,
  - zone-specific diffusivity modeling,
  - and coupled transport calculations.

---

# Physical System

The system consists of:
- a porous catalytic ammonia decomposer,
- a hydrogen–air combustor,
- and an annular flue-gas recirculation region used to thermally sustain the decomposition process.

The combustor exhaust is partially recirculated through the annular decomposer region to supply the heat required for the endothermic ammonia cracking reaction.

---

# Numerical Methodology

The CFD framework includes:
- 2D axisymmetric reacting-flow simulations,
- Brinkman–Darcy–Forchheimer porous media modeling,
- finite-rate Arrhenius chemistry,
- species transport,
- laminar diffusion flame modeling,
- conjugate heat transfer,
- and pressure-based steady-state formulation.

ANSYS Fluent was used for the CFD simulations, while MATLAB was used for reduced-order analysis, plug-flow modeling, and thermal recirculation calculations.

---

# Repository Structure

```text
Presentation/        -> Thesis presentation slides
MATLAB Analysis/     -> MATLAB scripts, validation, and post-processing
UDFs/                -> ANSYS Fluent User-Defined Functions
```

---

# MATLAB Scripts

### `Chein_ODE_MATLAB.m`
Derives and solves the plug-flow species balance governing ammonia decomposition by retaining the full convective flux term and recovering the initial value problem (IVP) formulation used by Chein et al.

The equations are solved numerically using an implicit backward-Euler discretization with nonlinear iteration at each axial step. The numerical solution closely reproduces the analytical implicit solution for NH₃, H₂, and N₂ mole fractions.

---

### `AmmoniaCombustor.m`
Processes and analyzes ANSYS Fluent combustor output data to extract thermodynamic, transport, and reaction-kinetics parameters required for the coupled decomposer–combustor framework.

The script computes quantities such as:
- inlet velocity,
- activation energy,
- pre-exponential factor,
- permeability,
- viscous resistance,
- and Forchheimer coefficient

used to reproduce and extend Chein’s 2D axisymmetric decomposer configuration.

The script also:
- solves the 1D plug-flow decomposition model using `fsolve`,
- validates extracted kinetic parameters,
- and evaluates porous-medium pressure drop using Ergun’s formulation.

---

### `Lean_CoupledAnalysis.m`
Performs thermal and energy-balance analysis of the thermally coupled decomposer–combustor system under lean combustion conditions (\(\lambda = 2.85\)).

The script computes:
- recirculation ratios,
- NTU-based heat-transfer parameters,
- flue-gas requirements,
- and decomposer operating conditions.

---

### `Stoichiometric_CoupledAnalysis.m`
Performs thermal and energy-balance analysis of the thermally coupled decomposer–combustor system under stoichiometric combustion conditions (\(\lambda = 1\)).

The script evaluates:
- recirculation requirements,
- NTU behavior,
- and thermal coupling

necessary to sustain ammonia decomposition near the target operating temperature.

---

### `PostProcessing.m`
Post-processing script used to compare the reproduced CFD results with digitized literature data from Chein et al.

The script overlays:
- species mole-fraction profiles,
- decomposition trends,
- and reactor-response behavior

to validate agreement with the literature model and confirm recovery of plug-flow behavior.

---

# UDF-Based Coupling Strategy

User-Defined Functions (UDFs) were implemented to:
- dynamically couple combustor outlet conditions to decomposer inlet conditions,
- compute recirculated flue-gas properties,
- and apply zone-specific transport models.

The UDF framework enables:
- closed-loop thermal interaction,
- iteration-wise recirculation updates,
- and stable coupled convergence.

---

# Zone-Specific Diffusivity Modeling

A key numerical challenge arose from conflicting transport requirements between the decomposer and combustor regions.

The decomposer required negligible axial diffusion to recover plug-flow behavior, while the combustor required physically realistic diffusion to sustain the hydrogen–air diffusion flame.

To resolve this:
- artificially low diffusivity was applied in the porous decomposer region,
- while Chapman–Enskog diffusivity was retained in the combustor region.

This was implemented using zone-specific UDFs.

The approach enabled simultaneous preservation of:
- plug-flow decomposition physics,
- and diffusion-controlled combustion physics

within the same computational domain.

---

# Key Results

- Successfully reproduced Chein et al.’s decomposition behavior
- Achieved near-complete NH₃ conversion
- Demonstrated stable thermal coupling
- Maintained decomposer operation near 893 K
- Achieved self-sustained coupled combustor–decomposer behavior
- Preserved both:
  - plug-flow decomposition physics
  - diffusion-controlled combustion physics

---

# Software Used

- ANSYS Fluent 2024 R1
- MATLAB

---
