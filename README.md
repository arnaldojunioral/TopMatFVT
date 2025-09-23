# TopMatFVT

This repository provides a MATLAB implementation of **topology optimization for periodic material microstructures**, based on the **Finite-Volume Theory (FVT)**. The algorithm computes the **homogenized constitutive matrix** of a periodic cell and optimizes its material distribution to achieve desired effective properties.

Supported interpolation models:
- **SIMP** (Solid Isotropic Material with Penalization)  
- **RAMP** (Rational Approximation of Material Properties)

Optional filtering techniques:
- **Sensitivity filter**  
- **Density filter**  

----

## 📌 Features

- Periodic boundary conditions automatically enforced.  
- Homogenization of effective elastic properties.  
- Objective function options:
  - **Shear modulus maximization**  
  - **Bulk modulus maximization**  
  - **Poisson’s ratio minimization**  
- Continuation scheme on penalization factors.  
- Initial material heterogeneity is defined by a circular void.  
- Optional **density/sensitivity filtering** for regularization.  

----

## ⚙️ Requirements

- MATLAB **R2015 or later**.
- GNU **Octave 6.4 or later**. 
- No extra toolboxes required.  

----

## 🚀 Getting started

Save the [TopMatFVT.m](https://raw.githubusercontent.com/arnaldojunioral/TopMatFVT/main/TopMatFVT.m) program (17 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

Run the main function:

**TopMatFVT(nx, ny, volfrac, penal, rfil, ft)**

where **nx** and **ny** define the number of subvolumes along the x- and y-directions, respectively; **volfrac** is the volume fraction constraint of solid material; **penal** is the penalization factor; and **rfil** and **ft** are additional parameters (filter radius and filter type) for the filtering analysis.

Run the main function:

```matlab

% In the following examples, the model consists of a structured mesh discretized into 100 × 100 subvolumes, with a solid material volume fraction constrained to 50%.

TopMatFVT(100, 100, 0.5, 3, [], []);      % Example with fixed penalization (penal factor = 3) and no filtering
TopMatFVT(100, 100, 0.5, 3, 2, 1);        % Example with fixed penalization (penal factor = 3) and sensitivity filter (rfil = 2 and ft = 1)
TopMatFVT(100, 100, 0.5, 1:3, 2, 2);      % Example with continuation scheme (penal factor = 1 to 3) and density filter (rfil = 2 and ft = 2)

````

The table below summarizes the key input parameters used in the simulation, including beam geometry, material properties, loading conditions, and visualization settings.

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

<!-- ## Documentation -->

<!-- The journal article uses the FVT3DELASTIC to generate the examples presented. -->

----

## 🎥 Topology evolution

<table align="center">
  <tr>
    <td align="center" valign="top">
      <strong>Shear modulus maximization</strong><br>
      <img width="300" height="300" alt="Shear modulus maximization" src="https://github.com/user-attachments/assets/2a03823e-f2d8-4459-9fa9-95ab7baacd67" />
    </td>
    <td align="center" valign="top">
      <strong>Bulk modulus maximization</strong><br>
      <img width="300" height="300" alt="Bulk modulus maximization" src="https://github.com/user-attachments/assets/869ee8cb-2cf9-4aba-896e-d710c8eded94" />
    </td>
    <td align="center" valign="top">
      <strong>Poisson ration minimization</strong><br>
      <img width="300" height="300" alt="Poisson ration minimization" src="https://github.com/user-attachments/assets/0756a229-cbe9-4301-9c0f-991c861e9062" />
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

The following table summarizes the four relevant references supporting the development of the proposed three-dimensional finite-volume theory. These works were selected based on their conceptual alignment with the present formulation, their methodological contributions, and their scientific impact.

| Rank | Reference                                                                                          | Relevance to the Study      | Scientific Impact | Justification                                                                                                                                                  |
|------|----------------------------------------------------------------------------------------------------|-----------------------------|-------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1    | Cavalcante, M.A.A., Pindera, M.-J. (2012a). *Generalized finite-volume theory for elastic analysis in solid mechanics.* Part I: Framework. *Journal of Applied Mechanics* | ⭐⭐⭐⭐⭐ | High              | Establishes the theoretical foundation of the generalized FVT used in this work. Introduces a novel framework that overcomes key challenges in solid mechanics modeling. |
| 2    | Cavalcante, M.A.A., Pindera, M.-J. (2012b). *Generalized finite-volume theory for elastic analysis in solid mechanics.* Part II: Results. *Journal of Applied Mechanics* | ⭐⭐⭐⭐⭐ | High              | Complements Part I by validating the FVT framework through numerical results, demonstrating its accuracy and robustness for linear elasticity problems.          |
| 3    | Cardiff, P., Demirdžić, I. (2021). *Thirty years of the finite volume method for solid mechanics.* *Archives of Computational Methods in Engineering* | ⭐⭐⭐⭐⭐ | Very High         | Offers a critical review of FVM developments, including FVT, positioning the current study within the broader trajectory of computational solid mechanics.       |
| 4    | Araujo, M.V.O., Lages, E.N., Cavalcante, M.A.A. (2020). *Checkerboard-free topology optimization for compliance minimization applying the finite-volume theory.* *Mechanics Research Communications* | ⭐⭐⭐⭐ | Medium            | Demonstrates the applicability of FVT beyond basic elasticity problems, showcasing its potential in advanced structural optimization scenarios.                |

