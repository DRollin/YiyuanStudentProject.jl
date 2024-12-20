```@meta
CurrentModule = YiyuanStudentProject
```

# Chemo-mechanical Problem

Documentation for [YiyuanStudentProject.jl](https://github.com/DRollin/YiyuanStudentProject.jl.git) 

The aim of this project is to develop and implement a linear transient chemo-mechanical multi-scale model within the FE² framework. To expedite the computation for a substantial number of Representative Volume Element (RVE) problems, the preparation for employing Numerical Model Reduction (NMR) utilizing snapshot Proper Orthogonal Decomposition (POD) shall be undertaken.

As indicated by the name, a chemo-mechanical problem refers to the situation where chemical reactions and mechanical changing are coupled together. Such problems are common in various fields, including battery science. Thus, a lithium-ion structural battery problem is discussed in this project, where the deformation and ion concentration have an influence of the chemical potential change.

This document provides an examination of the mathematical model and comprehensive explanation of the code implementation with an example.