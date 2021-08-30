# Contents

[Link 1](# Available-solvers-with-hybrid-approximation)

# Available solvers with hybrid approximation

United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions.

Only OpenFOAM+ version of the OpenFOAM technology is supported sinces 2018. Use branches digitef-dev-YYMM, where YYMM - corresponds to OpenFOAM+ version, for example 1812

1. **pimpleCentralFoam** - Pressure-based semi implicit compressible flow of perfect gas solver based on 
central-upwind schemes of Kurganov and Tadmor with LTS support for steady-state calculations

2. **rhoPimpleCentralFoam** - Pressure-based semi implicit compressible flow  of real gas solver based on 
central-upwind schemes of Kurganov and Tadmor and LTS support for steady-state calculations

3. **pimpleCentralDyMFoam** - Pressure-based semi implicit compressible flow of perfect gas solver based 
on central-upwind schemes of Kurganov and Tadmor with mesh motion and LTS support for steady-state calculations

4. **reactingPimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of 
Kurganov and Tadmor for combustion with chemical reactions and LTS support for steady-state calculations

5. **twoPhaseMixingCentralFoam** - Transient Eulerian two-phase solver. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

6. **twoPhaseMixingCentralDyMFoam** - Transient Eulerian two-phase solver with dynamic meshes. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

7. **chtMultiRegionCentralFoam** -     Pressure-based semi implicit solver, based on hybrid central-upwind schemes
    of Kurganov and Tadmor for conjugate simulation of compressible flows of perfect gas (Mach 
    number is ranging from 0 to 6) and solid body heat transfer.

# Discussion Telegram group

You can discuss questions of [hybridCentralSolvers](https://github.com/unicfdlab/hybridCentralSolvers) usage at Telegram Group: https://t.me/hybridCentralSolvers 

Available OpenFOAM versions:
* OpenFOAM 3.1 - master branch
* OpenFOAM 4.1 - dev-of4.1 branch
* OpenFOAM 6   - dev-of6 branch
* OpenFOAM+ 1812 - digitef-dev-1812
* OpenFOAM+ 1912 - digitef-dev-1912
* OpenFOAM+ 2012 - digitef-dev-2012

## The library/approach was useful in next research studies.

### >>>>> 2021 <<<<<

| Title | Description |
|------|-------------|
|[The Eulerian–Lagrangian Approach for the Numerical Investigation of an Acoustic Field Generated by a High-Speed Gas-Droplet Flow ](https://www.mdpi.com/2311-5521/6/8/274)| ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg) |
|[Dynamics of detonation transmission and propagation in a curved chamber: a numerical and experimental analysis](https://doi.org/10.1016/j.combustflame.2020.09.032)|![Experiment vs calculation](https://ars.els-cdn.com/content/image/1-s2.0-S0010218020304168-gr2.jpg)|
|[Modelling of Supersonic and Subsonic Flows Using Hybrid PressureBased Solver in Openfoam](https://doi.org/10.11159/ffhmt21.107)|![Schematic of Bluff body burner](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Bluff-body-sketch.png)|

### >>>>> 2020 <<<<<

| Title | Description |
|------|-------------|
|[Numerical simulation of transpiration cooling experiments in supersonic flow using OpenFOAM](https://link.springer.com/article/10.1007/s12567-019-00292-6)|![Schematic illustration of the applied porous interface model](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs12567-019-00292-6/MediaObjects/12567_2019_292_Fig1_HTML.png?as=webp)|
|[On uncertainty quantification via the ensemble of independent numerical solutions](https://doi.org/10.1016/j.jocs.2020.101114)|![Flow scheme](https://ars.els-cdn.com/content/image/1-s2.0-S1877750319310695-gr1.jpg)|
|[Influence of hydrogen equivalence ratios on supersonic combustion based on large eddy simulations](https://doi.org/10.1016/j.ijhydene.2020.02.054)|![Model scramjet](https://ars.els-cdn.com/content/image/1-s2.0-S036031992030584X-gr1.jpg)|
|[A quasi-direct numerical simulation solver for compressible reacting flows](https://doi.org/10.1016/j.compfluid.2020.104718)|![Data exchange between OpenFOAM and Cantera](https://ars.els-cdn.com/content/image/1-s2.0-S0045793020302887-gr1.jpg)|
|[ON THE CONSTRUCTION OF  A GENERALIZED COMPUTATIONAL EXPERIMENT IN VERIFICATION PROBLEMS](https://lppm3.ru/files/journal/XLVIII/MathMontXLVIII-Alekseev.pdf)|![#D flow around cone](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Cone-exact-3D.png)|
|[Enhanced Pressure Based Coupled Algorithm to Combine with Pressure–Velocity-Enthalpy for all Mach Number Flow](https://link.springer.com/article/10.1007/s42405-020-00337-9)|![Original pressure based p–h coupled algorithm](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs42405-020-00337-9/MediaObjects/42405_2020_337_Fig1_HTML.png?as=webp)|
|[Pressure-Based Solution Framework for Non-Ideal Flows at All Mach Numbers](https://link.springer.com/chapter/10.1007/978-3-030-49626-5_4)|![Fully developed jet structure of a n-hexane jet injected into a quiescent nitrogen atmosphere](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-49626-5_4/MediaObjects/492738_1_En_4_Fig3_HTML.png)|
|[Entwicklung eines Simulationsmodells für Schaltlichtbögen in Überspannungsableitern](https://publikationsserver.tu-braunschweig.de/servlets/MCRFileNodeServlet/dbbs_derivate_00047899/Diss_Sander_Christian.pdf)|![Spark](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/spark-tubraunscheig.png)|

###  >>>>> 2019 <<<<<
| Title | Description |
|------|-------------|
|[Numerical method to simulate detonative combustion of hydrogen-air mixture in a containment](https://doi.org/10.1080/19942060.2019.1660219)| ![Containement](https://www.tandfonline.com/na101/home/literatum/publisher/tandf/journals/content/tcfm20/2019/tcfm20.v013.i01/19942060.2019.1660219/20191106/images/medium/tcfm_a_1660219_f0018_oc.jpg)|
|[Numerical investigation of the auto-ignition of transient hydrogen injection in supersonic airflow](https://doi.org/10.1016/j.ijhydene.2019.07.215)|![Shadowgraph of the jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/1-s2.0-S0360319919328514-gr4.jpg)|
|[Verification on the Ensemble of Independent Numerical Solutions](https://link.springer.com/chapter/10.1007/978-3-030-22750-0_25)|![Solvers comparison](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-22750-0_25/MediaObjects/485772_1_En_25_Fig2_HTML.png)|
|[Computational Study of Reactants Mixing in a Rotating Detonation Combustor Using Compressible RANS](https://link.springer.com/article/10.1007/s10494-019-00097-x)|![Detailed shock structure of the baseline flow case in the injection region derived from the Mach number contour plot at the longitudinal mid-plane](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs10494-019-00097-x/MediaObjects/10494_2019_97_Fig5_HTML.png?as=webp)|

### >>>>> 2018 <<<<<
| Title | Description |
|------|-------------|
|[Real-Gas Effects and Phase Separation in Underexpanded Jets at Engine-Relevant Conditions](https://doi.org/10.2514/6.2018-1815)|![Jet development history](https://www.researchgate.net/profile/Christoph-Traxinger/publication/322309300/figure/fig5/AS:622107033612289@1525333283581/figure-fig5_W640.jpg)|
|[Analysis of the Accuracy of OpenFOAM Solvers for the Problem of Supersonic Flow Around a Cone](https://link.springer.com/chapter/10.1007/978-3-319-93713-7_18)|![Cone sketch](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-319-93713-7_18/MediaObjects/469704_1_En_18_Fig1_HTML.gif)|
|[Development of a new OpenFOAM solver using regularized gas dynamic equations](https://doi.org/10.1016/j.compfluid.2018.02.010)|![Ladenburgh jet](https://ars.els-cdn.com/content/image/1-s2.0-S0045793018300641-gr14.jpg)|
|[A hybrid pressure-based solver for nonideal single-phase fluid flows at all speeds](https://doi.org/10.1002/fld.4512)|![Experiment vs. calculation](https://onlinelibrary.wiley.com/cms/asset/16108f70-6fec-4197-9aab-f84cbc5c2a1d/fld4512-fig-0005-m.jpg)|
|[Comparison of the Performance of Open-Source and Commercial CFD Packages for Simulating Supersonic Compressible Jet Flows](https://doi.org/10.1109/IVMEM.2018.00019)|![Airbag computational domain sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Airbag.png)|
|[Numerical modelling of two-dimensional perfect gas flows using RKDG method on unstructured meshes](https://doi.org/10.1063/1.5065323)|![RKDG (a) vs. rhoPimpleCentralFoam (b)](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/forwardStep-RKDG-vs-RPCF.png)|
|[CFD methodologies for compressible atomising and cavitating multi-phase flows](https://eprints.utas.edu.au/28677/)|![Jet iso contours](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Jet-Hongjiang.png)|

### >>>>> 2017 <<<<<
| Title | Description |
|------|-------------|
|[Numerical investigation on an array of Helmholtz resonators for the reduction of micro-pressure waves in modern and future high-speed rail tunnel systems](https://doi.org/10.1016/j.jsv.2017.04.022)| ![Helmholtz resonantors array mesh](https://ars.els-cdn.com/content/image/1-s2.0-S0022460X17303280-gr10.jpg)|
|[Comparative Study of the Accuracy for OpenFOAM Solvers](https://doi.org/10.1109/ISPRAS.2017.00028)|![Flow around the cone](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/index.jpeg)|
|[ Analysis of Radiation Discretization for Modelling a Spark Gap for Surge Currents ](https://doi.org/10.14311/ppt.2017.1.56)|![Spark sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Plasma-spark-sketch.png)|
|[Numerical analysis of cavitation about marine propellers using a compressible multiphase VOF fractional step method](https://www.researchgate.net/publication/319306852_Numerical_analysis_of_cavitation_about_marine_propellers_using_a_compressible_multiphase_VOF_fractional_step_method)|![Comparison of cavitation morphology between cFSMVOF and other commercial/open-source codes](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/cFSMVOF.png)|
|[Computational analysis and mitigation of micro-pressure waves in high-speed train tunnels](https://doi.org/10.25560/72653)|![Typical high speed train](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Train-1.png)|
|[Numerical study of characteristic modes and frequencies of flow in high-speed compressors (in English)](https://doi.org/10.15514/ISPRAS-2017-29(1)-2)|![ERCOFTAC Pump sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/ERCOFTAC-pump.png)|

###  >>>>> 2016 <<<<<
| Title | Description |
|------|-------------|
|[On the Stability of Supersonic Boundary Layers with Injection](https://thesis.library.caltech.edu/9755/)| ![Scheme of boundary layer interaction with jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/boundary-layer.png)|
|[Study of capabilities of hybrid scheme for advection terms approximation in mathematical models of compressible flows (in Russian)](https://ispranproceedings.elpub.ru/jour/article/view/121)|![Liquid ring vacuum pump](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/LRVP.png)|
|[LES Discretization Methods for Unstructured Meshes Based on the Finite Volume Method](https://doi.org/10.21656/1000-0887.370228)|![Vorticity: Flow around cylinder](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Vorticity-vs-scheme.png)|

   When using these solvers, please cite the following works:
   * [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878441.svg)](https://doi.org/10.5281/zenodo.3878441)
   * Kraposhin MV, Banholzer M, Pfitzner M, Marchevsky IK. A hybrid pressure‐based solver for nonideal single‐phase fluid flows at all speeds. Int J Numer Meth Fluids. 2018;88:79–99. https://doi.org/10.1002/fld.4512
   * Kraposhin MV, Strijhak SV, Bovtrikova A Adaptation of Kurganov-Tadmor Numerical Scheme for Applying in Combination with the PISO Method in Numerical Simulation of Flows in a Wide Range of Mach Numbers. Procedia Computer Science. 2015;66:43-52. https://doi.org/10.1016/j.procs.2015.11.007
