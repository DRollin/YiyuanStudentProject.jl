```@meta
CurrentModule = YiyuanStudentProject
```

# Chemo-Mechanical Problem

Documentation for [YiyuanStudentProject.jl](https://github.com/DRollin/YiyuanStudentProject.jl.git) 

## Introduction

The aim of this project is to develop and implement a linear transient chemo-mechanical multi-scale model within the FE²-framework. Representitive Volumen Element (RVE) with fully resolved material microstructure is established on the sub scale, whereby a homogeneous macro scale is assumed and governed by the corresponding sub-scale behaviour.

As indicated by the name, a chemo-mechanical problem refers to the situation where transport of chemical species and deformation are coupled. Such problems are common in various fields, including battery science. Thus, a lithium-ion structural battery problem is discussed in this project, where the deformation and ion concentration have an influence on the chemical potential.

This document provides an investigation of the mathematical model and comprehensive explanation of the implementation with an example case study. 

## Declaration of AI use

The following AI tools are used helping analysing code error messages and improving structure of the documentation as well as proofreading the documentation:
- ChatGPT 4o
- Microsoft Copilot
