# ✨ TopMatFVT

This repository provides a MATLAB implementation of **topology optimization for periodic material microstructures**, based on the **Finite-Volume Theory (FVT)**. The algorithm computes the **homogenized constitutive matrix** of a periodic cell and optimizes its material distribution to achieve desired effective properties.

Supported material interpolation models:
- **SIMP** (Solid Isotropic Material with Penalization)  
- **RAMP** (Rational Approximation of Material Properties)

Optional filtering techniques:
- **Sensitivity filter**  
- **Density filter**  

----

## 📌 Features

- Periodic boundary conditions automatically enforced.  
- Homogenization of periodic cellular materials based on the concept of Repeating Unit Cell (RUC).  
- Objective function options:
  - **Shear modulus maximization**  
  - **Bulk modulus maximization**  
  - **Poisson’s ratio minimization**  
- Either a continuation scheme applied to penalization factors or a fixed penalization approach.
- Initial material heterogeneity is defined by a circular void.  
- Optional **sensitivity/density filtering** for the solution's regularization.  

----

## ⚙️ Requirements

The implementation is fully compatible with both **MATLAB (R2015 or later)** and **GNU Octave (version 6.4 or later)**. No additional toolboxes are required, ensuring that the code can be executed in a standard installation of either environment.

----

## 🚀 Getting started

Save the [TopMatFVT.m](https://raw.githubusercontent.com/arnaldojunioral/TopMatFVT/main/TopMatFVT.m) program (9.57 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

**TopMatFVT(nx, ny, volfrac, penal, rfil, ft)**

where **nx** and **ny** define the number of subvolumes along the x- and y-directions, respectively; **volfrac** is the volume fraction constraint of solid material; **penal** is the penalization factor (fixed or continuation scheme); **rfil** and **ft** are additional parameters (filter radius and filter type) for the filtering analysis.

Run the main function:

```matlab

% In the following examples, the model consists of a structured mesh discretized into 100 × 100 subvolumes, with a solid material volume fraction constrained to 50%.

TopMatFVT(100, 100, 0.5, 3, [], []);      % Example with fixed penalization (penal factor = 3) and no filtering solution (..., [], [])
TopMatFVT(100, 100, 0.5, 3, 2, 1);        % Example with fixed penalization (penal factor = 3) and sensitivity filter (filter radius rfil = 2 and filter type ft = 1)
TopMatFVT(100, 100, 0.5, 1:3, 2, 2);      % Example with continuation scheme (penal factor = 1 to 3) and density filter (filter radius rfil = 2 and filter type ft = 2)

````

The table below summarizes the main input parameters considered in the simulations, including the material properties, model settings, objective function, and numerical update controls.

##### Model Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| E0        | 1.0   | Young's modulus of solid material |
| nu        | 0.3   | Poisson's ratio |
| ctp       | 3     | Objective function: 1 (shear modulus), 2 (bulk modulus), 3 (Poisson's ratio) |
| R         | min(nx,ny)/6 | Radius of circular material heterogeneity |
| mdl       | 'SIMP' | Material interpolation method: 'SIMP' or 'RAMP' |
| eta       | 1/3   | Damping factor |
| move      | 0.2   | Move limit for design variable update |

## 📘 Documentation

For further details on the theoretical background and verification, please refer to the following articles:

- **Santos Júnior, A.** & **Cavalcante, M. A. A.**, (2025). *Checkerboard-free topology optimization for cellular materials via the finite-volume theory*. *Optimization and Engineering*. [https://doi.org/10.1007/s11081-025-09988-7](https://doi.org/10.1007/s11081-025-09988-7)

- **Santos Júnior, A.** & **Cavalcante, M. A. A.**, (2024). *Topology optimization of periodic cellular materials employing the finite-volume theory*. *Engineering Optimization*. [https://doi.org/10.1080/0305215X.2024.2379019](https://doi.org/10.1080/0305215X.2024.2379019)

----

## 🎥 Topology evolution

<p align="center"><strong> No filtering solutions</strong></p>

<table align="center">
  <tr>
    <td align="center" valign="top">
      <strong>Shear modulus maximization</strong><br>
      <img width="250" height="250" alt="Shear modulus maximization" src="https://github.com/user-attachments/assets/2a03823e-f2d8-4459-9fa9-95ab7baacd67" />
    </td>
    <td align="center" valign="top">
      <strong>Bulk modulus maximization</strong><br>
      <img width="250" height="250" alt="Bulk modulus maximization" src="https://github.com/user-attachments/assets/869ee8cb-2cf9-4aba-896e-d710c8eded94" />
    </td>
    <td align="center" valign="top">
      <strong>Poisson ratio minimization</strong><br>
      <img width="250" height="250" alt="Poisson's ratio minimization" src="https://github.com/user-attachments/assets/0756a229-cbe9-4301-9c0f-991c861e9062" />
    </td>
  </tr>
</table>

----

## Error reporting

We strive to ensure that the implementation of the finite-volume theory is accurate, efficient, and well-documented. However, if you encounter unexpected behavior, inconsistencies, or potential bugs in the code, we welcome your feedback.

Please feel free to:

- **Open an issue** on the [GitHub repository](https://github.com/arnaldojunioral/FVT3DELASTIC/issues) describing the problem in detail.
- Or **contact us directly** via email at [arnaldo@ctec.ufal.br](mailto:arnaldo@ctec.ufal.br).

Your contributions help improve the reliability and usability of this project for the research community.

----

## Authors

Project developed by:

* Arnaldo dos Santos Júnior  arnaldo@ctec.ufal.br
* Márcio André Araújo Cavalcante marcio.cavalcante@ceca.ufal.br

----

## References

The following table summarizes the six relevant references supporting the development of the proposed three-dimensional finite-volume theory. These works were selected based on their conceptual alignment with the present formulation, their methodological contributions, and their scientific impact.

| Rank | Reference | Relevance to the Study | Scientific Impact | Justification |
|------|-----------|----------------------|-----------------|---------------|
| 1    | Cavalcante, M.A.A., Pindera, M.-J. (2012a). *Generalized finite-volume theory for elastic stress analysis in solid mechanics.* Part I: Framework. *Journal of Applied Mechanics, Transactions ASME*, v. 79, p. 051006.| ⭐⭐⭐⭐⭐ | High | Establishes the theoretical foundation of the generalized FVT used in this work, introducing a novel framework that overcomes key challenges in solid mechanics modeling. |
| 2    | Cavalcante, M.A.A., Pindera, M.-J.; Khatam, H. (2012). *Finite-volume micromechanics of periodic materials: Past, present and future.* *Composites: Part B*, v. 43, p. 2521–2543. | ⭐⭐⭐⭐⭐ | High | Reviews FVT for periodic materials, summarizing theory, applications, and state-of-the-art, directly relevant for model validation. |
| 3    | Bendsøe, M. P.; Kikuchi, N. (1988). *Generating optimal topologies in structural design using a homogenization method.* *Computer Methods in Applied Mechanics and Engineering*, v. 71, p. 197–224. | ⭐⭐⭐⭐⭐ | Very High | Classic work in topology optimization using homogenization, introducing the approach foundational to modern topology optimization. |
| 4    | Bendsøe, M. P.; Sigmund, O. (1999). *Material interpolation schemes in topology optimization.* *Archive of Applied Mechanics*, v. 69, p. 635–654. | ⭐⭐⭐⭐ | High | Discusses SIMP and general material interpolation schemes central to penalization strategies in topology optimization. |
| 5    | Santos Júnior, A.; Cavalcante, M. A. A. (2025). *Checkerboard-free topology optimization for cellular materials via the finite-volume theory.* *Optimization and Engineering*. | ⭐⭐⭐ | Medium | Provides a comparative study between energy equivalence and mean-field theory approaches to topology optimization of periodic cellular materials. |
| 6    | Santos Júnior, A.; Cavalcante, M. A. A. (2024). *Topology optimization of periodic cellular materials employing the finite-volume theory.* *Engineering Optimization*. | ⭐⭐⭐ | Medium | Extends FVT-based topology optimization to periodic cellular microstructures employing the energy equivalence. |
| 7    | Araujo, M. V. O.; Lages, E. N.; Cavalcante, M. A. A. (2020a). *Checkerboard-free topology optimization for compliance minimization applying the finite-volume theory.* *Mechanics Research Communications*, v. 108, p. 103581. | ⭐⭐⭐ | Medium | Demonstrates checkerboard-free topology optimization using FVT for compliance minimization, complementing modern applications. |


