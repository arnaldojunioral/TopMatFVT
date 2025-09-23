# TopMatFVT

This repository provides a **free** MATLAB implementation of a Topology Optimization of Periodic Materials using the Standard Finite-Volume Theory (FVT) Formulation.. The code offers:

* an efficient and flexible code for conducting numerical investigations of periodic cellular materials;
* stress and displacement analyses through a local equilibrium-based formulation; and
* easy extension to other problems involving three-dimensional structural components.

## Getting started

Save the [TopMatFVT.m](https://raw.githubusercontent.com/arnaldojunioral/TopMatFVT/main/TopMatFVT.m) program (17 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

**TopMatFVT(n1, n2, n3)**

where **n1**, **n2**, and **n3** define the number of subvolumes along the x<sub>1</sub>, x<sub>2</sub>, and x<sub>3</sub> directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below.
<!-- <p align="center">
<img width="350" height="350" alt="image" src="https://github.com/user-attachments/assets/3d92838e-2fcb-40f7-b0da-80d891ec62d6" />
</p> -->
<p align="center">
<img width="510" height="280" alt="image" src="https://github.com/user-attachments/assets/9efeaf36-5ae0-45d4-b3fc-7a83c7a8b952" />
</p>

The table below summarizes the key input parameters used in the simulation, including beam geometry, material properties, loading conditions, and visualization settings.

#### Model Parameters

| Parameter       | Description                                          | Unit             |
|-----------------|------------------------------------------------------|------------------|
| L             | Beam length                                          | mm               |
| H             | Beam height                                          | mm               |
| B             | Beam width                                           | mm               |
| E             | Young's modulus (material stiffness)                 | MPa              |
| nu            | Poisson's ratio                                      | –                |
| P             | Applied load (negative indicates downward force)     | N                |
| pb            | Problem: 'flexure', 'torsion', or 'torsion-flexure' | –          |
| amp           | Amplification factor for deformation visualization   | –                | 

<!-- ## Documentation -->

<!-- The journal article uses the FVT3DELASTIC to generate the examples presented. -->

## Example: 3D Cantilever beam

This example presents a standard benchmark problem involving a three-dimensional cantilever beam. The analysis domain and boundary conditions, shown in the figure below, are used to verify the functionality of the FVT3DELASTIC code.

<p align="center">
<img width="400" height="206" alt="image" src="https://github.com/user-attachments/assets/56962135-65c7-44d0-ba56-1cc5c18a9910" />
</p>

The parameters used in the analysis are listed in the table below. The function call **FVT3DELASTIC(75, 16, 16)** corresponds to a structured mesh discretized into 75 × 16 × 16 subvolumes.

| L   | H   | B   | E      | ν   | P    | pb       | amp   |
|-----|-----|-----|--------|-----|------|----------|-------|
| 500 | 100 | 100 | 150000 | 0.3 | 2000 | 'flexure'  | 1000  |

#### Output

<table align="center">
  <tr>
    <td align="center" valign="top">
      <strong>Deformed structure</strong><br>
      <img width="370" height="200" alt="Deformed" src="https://github.com/user-attachments/assets/55c63898-61b3-42af-91ae-f9514570bcb1" />
    </td>
    <td align="center" valign="top">
      <strong>Stress fields</strong><br>
      <img width="900" height="300" alt="Stress fields" src="https://github.com/user-attachments/assets/1a826ba6-0fc4-4749-b6d4-4bfef5d56534" />
    </td>
  </tr>
</table>

## Error reporting

We strive to ensure that the implementation of the finite-volume theory is accurate, efficient, and well-documented. However, if you encounter unexpected behavior, inconsistencies, or potential bugs in the code, we welcome your feedback.

Please feel free to:

- **Open an issue** on the [GitHub repository](https://github.com/arnaldojunioral/FVT3DELASTIC/issues) describing the problem in detail.
- Or **contact us directly** via email at [arnaldo@ctec.ufal.br](mailto:arnaldo@ctec.ufal.br).

Your contributions help improve the reliability and usability of this project for the research community.

## Authors

Project developed by:

* Arnaldo dos Santos Júnior  arnaldo@ctec.ufal.br
* Márcio André Araújo Cavalcante marcio.cavalcante@ceca.ufal.br

## References

The following table summarizes the four relevant references supporting the development of the proposed three-dimensional finite-volume theory. These works were selected based on their conceptual alignment with the present formulation, their methodological contributions, and their scientific impact.

| Rank | Reference                                                                                          | Relevance to the Study      | Scientific Impact | Justification                                                                                                                                                  |
|------|----------------------------------------------------------------------------------------------------|-----------------------------|-------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1    | Cavalcante, M.A.A., Pindera, M.-J. (2012a). *Generalized finite-volume theory for elastic analysis in solid mechanics.* Part I: Framework. *Journal of Applied Mechanics* | ⭐⭐⭐⭐⭐ | High              | Establishes the theoretical foundation of the generalized FVT used in this work. Introduces a novel framework that overcomes key challenges in solid mechanics modeling. |
| 2    | Cavalcante, M.A.A., Pindera, M.-J. (2012b). *Generalized finite-volume theory for elastic analysis in solid mechanics.* Part II: Results. *Journal of Applied Mechanics* | ⭐⭐⭐⭐⭐ | High              | Complements Part I by validating the FVT framework through numerical results, demonstrating its accuracy and robustness for linear elasticity problems.          |
| 3    | Cardiff, P., Demirdžić, I. (2021). *Thirty years of the finite volume method for solid mechanics.* *Archives of Computational Methods in Engineering* | ⭐⭐⭐⭐⭐ | Very High         | Offers a critical review of FVM developments, including FVT, positioning the current study within the broader trajectory of computational solid mechanics.       |
| 4    | Araujo, M.V.O., Lages, E.N., Cavalcante, M.A.A. (2020). *Checkerboard-free topology optimization for compliance minimization applying the finite-volume theory.* *Mechanics Research Communications* | ⭐⭐⭐⭐ | Medium            | Demonstrates the applicability of FVT beyond basic elasticity problems, showcasing its potential in advanced structural optimization scenarios.                |

