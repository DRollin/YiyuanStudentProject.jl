var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference overview","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference overview","title":"Reference overview","text":"Pages = [\n    \"types.md\",\n    \"fine_scale.md\",\n]","category":"page"},{"location":"math_model/fine_scale/#Fine-Scale-Model","page":"Fine Scale Model","title":"Fine Scale Model","text":"","category":"section"},{"location":"math_model/fine_scale/#Strong-From","page":"Fine Scale Model","title":"Strong From","text":"","category":"section"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"Three primary variable fields: displacement vector (u(xt) Omega times mathbbR^+ rightarrow mathbbR^3) and molar ion concentration (c(xt) Omega times mathbbR^+ rightarrow mathbbR) as well as the chemical potential gradient (mu(xt) Omega times mathbbR^+ rightarrow mathbbR) are solved in the lithium structural battery problem. To solve the problem three equilibrium conditions in mechanics, diffusion and chemical potential are given in Strong Forms as following:","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"-sigma cdot nabla = 0\n\nin  Omega times (0T \n\ndotc + j cdot nabla = 0\nin Omega times (0T\n\nmu - mu^en = 0\nin Omega times (0T ","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"Where With the boundary conditions, where the boundary is divied into three parts for the three variable fields representivly  (Gamma = Gamma_D^(u) cup Gamma_N^(u) = Gamma_D^(mu) cup Gamma_N^(mu) = Gamma_D^(c) cup Gamma_N^(c)):","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"u = u^p\non Gamma_D^(u) times (0T\n\nt = sigma cdot n = t^p\non Gamma_N^(u) times (0T\n\nmu = mu^p\non Gamma_D^(mu) times (0T \n\nh = -j cdot n = h^p\non Gamma_N^(mu) times (0T","category":"page"},{"location":"math_model/fine_scale/#Constitutive-Models","page":"Fine Scale Model","title":"Constitutive Models","text":"","category":"section"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"Due to the coupling of the mechanical and chemical aspects the total strain has both the contribution due to a deformation state u(x,t) and the ion concentration c(x,t), where an ion intercalation tensor (alpha^ch) is used based on the models of xxxx and xxxxx. ","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"    epsilonu = (u otimes nabla)^sym\n\n    epsilon^ch(c) = alpha^ch (c-c_ref)","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"The free energy psi is assumed to be the sum of mechanical and chemical parts psi (epsilon c) = psi^mech(epsilonc) + psi^chem(c).","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"\n    psi^mech = frac12(epsilon-epsilon^ch(c)) colon E colon (epsilon-epsilon^ch(c))\n\n    psi^chem = (c-c_ref)mu_ref + frac12 fracRtheta_refc_m(c-c_ref)^2\n","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"Furthermore, a constant mobility tensor M for the assumption of a linear relation between ion flux and the gradient of the chemical potential is used using a mobility coefficient eta .","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"In this project a reference temperature theta_ref and the converntration c_ref as well as the reference converntration c_ref are gloable constant material parameters for the purpose of linearization. So that the simplified constitutive equations are:","category":"page"},{"location":"math_model/fine_scale/","page":"Fine Scale Model","title":"Fine Scale Model","text":"\n    sigma(epsilonc) = fracpartialpsipartialepsilon = E colon (epsilon-epsilon^ch(c))\n\n    mu^en(epsilonc) = fracpartialpsipartialepsilon = mu_ref + fracRtheta_refc_m(c-c_ref) - alpha^ch colon sigma(epsilonc)\n\n    j(nablamu) = -M cdot nablamu\n","category":"page"},{"location":"reference/fine_scale/","page":"Fine scale","title":"Fine scale","text":"CurrentModule = YiyuanStudentProject","category":"page"},{"location":"reference/fine_scale/#Fine-scale","page":"Fine scale","title":"Fine scale","text":"","category":"section"},{"location":"reference/fine_scale/","page":"Fine scale","title":"Fine scale","text":"generate_rve_grid\nget_volume\nget_phase_bias\nprepare_setup\nadd_bc!\nassemble_K!\nassemble_Kₑ!\nassemble_M!\nassemble_Mₑ!\ncompute_time_step!\nsolve_load_case\nplot_grid\nanimate_result","category":"page"},{"location":"reference/fine_scale/#YiyuanStudentProject.generate_rve_grid","page":"Fine scale","title":"YiyuanStudentProject.generate_rve_grid","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.get_volume","page":"Fine scale","title":"YiyuanStudentProject.get_volume","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.get_phase_bias","page":"Fine scale","title":"YiyuanStudentProject.get_phase_bias","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.prepare_setup","page":"Fine scale","title":"YiyuanStudentProject.prepare_setup","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.add_bc!","page":"Fine scale","title":"YiyuanStudentProject.add_bc!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.assemble_K!","page":"Fine scale","title":"YiyuanStudentProject.assemble_K!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.assemble_Kₑ!","page":"Fine scale","title":"YiyuanStudentProject.assemble_Kₑ!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.assemble_M!","page":"Fine scale","title":"YiyuanStudentProject.assemble_M!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.assemble_Mₑ!","page":"Fine scale","title":"YiyuanStudentProject.assemble_Mₑ!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.compute_time_step!","page":"Fine scale","title":"YiyuanStudentProject.compute_time_step!","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.solve_load_case","page":"Fine scale","title":"YiyuanStudentProject.solve_load_case","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.plot_grid","page":"Fine scale","title":"YiyuanStudentProject.plot_grid","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"reference/fine_scale/#YiyuanStudentProject.animate_result","page":"Fine scale","title":"YiyuanStudentProject.animate_result","text":"TODO\n\n\n\n\n\n","category":"function"},{"location":"math_model/#Mathmatical-Model","page":"Mathmatical Model overview","title":"Mathmatical Model","text":"","category":"section"},{"location":"math_model/","page":"Mathmatical Model overview","title":"Mathmatical Model overview","text":"Pages = [\n    \"fine_scale.md\",\n    \"upscaling.md\",\n    \"macro_scale.md\",\n]","category":"page"},{"location":"reference/types/","page":"Types","title":"Types","text":"CurrentModule = YiyuanStudentProject","category":"page"},{"location":"reference/types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"reference/types/","page":"Types","title":"Types","text":"Material\nRVE\nLoadCase\nPhaseSetup\nRVESetup","category":"page"},{"location":"reference/types/#YiyuanStudentProject.Material","page":"Types","title":"YiyuanStudentProject.Material","text":"TODO\n\n\n\n\n\n","category":"type"},{"location":"reference/types/#YiyuanStudentProject.RVE","page":"Types","title":"YiyuanStudentProject.RVE","text":"TODO\n\n\n\n\n\n","category":"type"},{"location":"reference/types/#YiyuanStudentProject.LoadCase","page":"Types","title":"YiyuanStudentProject.LoadCase","text":"TODO\n\n\n\n\n\n","category":"type"},{"location":"reference/types/#YiyuanStudentProject.PhaseSetup","page":"Types","title":"YiyuanStudentProject.PhaseSetup","text":"TODO\n\n\n\n\n\n","category":"type"},{"location":"reference/types/#YiyuanStudentProject.RVESetup","page":"Types","title":"YiyuanStudentProject.RVESetup","text":"TODO\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = YiyuanStudentProject","category":"page"},{"location":"#Chemo-mechanical-Problem","page":"Home","title":"Chemo-mechanical Problem","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for YiyuanStudentProject.jl ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The aim of this project is to develop and implement a linear transient chemo-mechanical multi-scale model within the FE² framework. To expedite the computation for a substantial number of Representative Volume Element (RVE) problems, the groundwork for employing Numerical Model Reduction (NMR) utilizing snapshot Proper Orthogonal Decomposition (POD) shall be undertaken.","category":"page"},{"location":"","page":"Home","title":"Home","text":"As can be told from the name, a chemo-mechanical problem refers to the situation where chemical reactions and mechanical changing are coupled together. Such problems are common in various fields, including battery science. Thus, a lithium-ion structural battery problem is discussed in this project. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This document provides an exhaustive examination of the mathematical model and comprehensive explanation the code implementation.","category":"page"}]
}
