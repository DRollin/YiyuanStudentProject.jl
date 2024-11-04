# Fine Scale Model

## Strong From

Three primary variable fields: displacement vector ($u(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}^{3}$) and molar ion concentration ($c(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$) as well as the chemical potential gradient ($\mu(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$) are solved in the lithium structural battery problem. To solve the problem three equilibrium conditions in mechanics, diffusion and chemical potential are given in Strong Forms as following:

```math
-\sigma \cdot \nabla = 0
\
in  \Omega \times (0,T] 
\\
\dot{c} + j \cdot \nabla = 0
in \Omega \times (0,T]
\\
\mu - \mu^{en} = 0
in \Omega \times (0,T] 
```

Where With the boundary conditions, where the boundary is divied into three parts for the three variable fields representivly  ($\Gamma = \Gamma_D^{(u)} \cup \Gamma_N^{(u)} = \Gamma_D^{(\mu)} \cup \Gamma_N^{(\mu)} = \Gamma_D^{(c)} \cup \Gamma_N^{(c)}$):

```math
u = u^{p}
on \Gamma_{D}^{(u)} \times (0,T]
\\
t := \sigma \cdot n = t^{p}
on \Gamma_{N}^{(u)} \times (0,T]
\\
\mu = \mu^{p}
on \Gamma_{D}^{(\mu)} \times (0,T] 
\\
h := -j \cdot n = h^{p}
on \Gamma_{N}^{(\mu)} \times (0,T]
```


## Constitutive Models
Due to the coupling of the mechanical and chemical aspects the total strain has both the contribution due to a deformation state u(x,t) and the ion concentration c(x,t), where an ion intercalation tensor ($\alpha^{ch}$) is used based on the models of xxxx and xxxxx. 

```math
    \epsilon[u] := (u \otimes \nabla)^{sym}
\\
    \epsilon^{ch}(c) := \alpha^{ch} (c-c_{ref})
```

The free energy $\psi$ is assumed to be the sum of mechanical and chemical parts $\psi (\epsilon, c) = \psi^{mech}(\epsilon,c) + \psi^{chem}(c)$.

```math

    \psi^{mech} := \frac{1}{2}(\epsilon-\epsilon^{ch}(c)) \colon E \colon (\epsilon-\epsilon^{ch}(c))
\\
    \psi^{chem} := (c-c_{ref})\mu_{ref} + \frac{1}{2} \frac{R\theta_{ref}}{c_m}(c-c_{ref})^2

```

Furthermore, a constant mobility tensor M for the assumption of a linear relation between ion flux and the gradient of the chemical potential is used using a mobility coefficient $\eta$ .

In this project a reference temperature $\theta_{ref}$ and the converntration $c_{ref}$ as well as the reference converntration $c_{ref}$ are gloable constant material parameters for the purpose of linearization. So that the simplified constitutive equations are:

```math

    \sigma(\epsilon,c) = \frac{\partial\psi}{\partial\epsilon} = E \colon (\epsilon-\epsilon^{ch}(c))
\\
    \mu^{en}(\epsilon,c) = \frac{\partial\psi}{\partial\epsilon} = \mu_{ref} + \frac{R\theta_{ref}}{c_m}(c-c_{ref}) - \alpha^{ch} \colon \sigma(\epsilon,c)
\\
    j(\nabla[\mu]) := -M \cdot \nabla[\mu]

```