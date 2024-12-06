# [Chemo-Mechanical Problem - fine scale](@id documentation-fine-scale)

## Strong From

Three primary variable fields: displacement vector ($u(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}^{3}$) 
and molar ion concentration ($c(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$) 
as well as the chemical potential gradient ($\mu(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$) 
are solved in the lithium structural battery problem. 
To solve the problem three equilibrium conditions in mechanics, diffusion and chemical potential are given in Strong Forms as following:

```math
\begin{align}
    -\sigma \cdot \nabla &= 0 & in  \Omega \times (0,T] 
\\
    \dot{c} + j \cdot \nabla &= 0 & in \Omega \times (0,T]
\\
    \mu - \mu^{en} &= 0 & in \Omega \times (0,T] 
\end{align}
```

Where the boundary is divied into three parts for the three variable fields representivly  
(``\Gamma = \Gamma_D^{(u)} \cup \Gamma_N^{(u)} = \Gamma_D^{(\mu)} \cup \Gamma_N^{(\mu)} = \Gamma_D^{(c)} \cup \Gamma_N^{(c)}``):

```math
\begin{align}
    u &= u^{p}
    &
    on \Gamma_{D}^{(u)} \times (0,T]
\\
    t := \sigma \cdot n &= t^{p}
    &
    on \Gamma_{N}^{(u)} \times (0,T]
\\
    \mu &= \mu^{p}
    &
    on \Gamma_{D}^{(\mu)} \times (0,T] 
\\
    h := -j \cdot n &= h^{p}
    &
    on \Gamma_{N}^{(\mu)} \times (0,T]
\end{align}
```

## Constitutive Models

Due to the coupling of the mechanical and chemical aspects the total strain has both the contribution due to a deformation state u(x,t)
and the ion concentration c(x,t), where an ion intercalation tensor ($\alpha^{ch}$) is used based on the models of xxxx and xxxxx. 

```math
\begin{align}
\epsilon[u] &:= (u \otimes \nabla)^{sym}
\\
\epsilon^{ch}(c) &:= \alpha^{ch} (c-c_{ref})
\end{align}
```

The free energy ``\psi`` is assumed to be the sum of mechanical and chemical parts 
``\psi (\epsilon, c) = \psi^{mech}(\epsilon,c) + \psi^{chem}(c)``.

```math
\begin{align}
    \psi^{mech} &:= \frac{1}{2}(\epsilon-\epsilon^{ch}(c)) \colon E \colon (\epsilon-\epsilon^{ch}(c))
\\
    \psi^{chem} &:= (c-c_{ref})\mu_{ref} + \frac{1}{2} \frac{R\theta_{ref}}{c_m}(c-c_{ref})^2
\end{align}
```

Furthermore, a constant mobility tensor M for the assumption of a linear relation between ion flux and the gradient of 
the chemical potential is used using a mobility coefficient $\eta$ .
In this project a reference temperature $\theta_{ref}$ and the converntration $c_{ref}$ as well as the
reference converntration $c_{ref}$ are gloable constant material parameters for the purpose of linearization. 
So that the simplified constitutive equations are:

```math
\begin{align}
    \sigma(\epsilon,c) &= \frac{\partial\psi}{\partial\epsilon}
    &= E \colon (\epsilon-\epsilon^{ch}(c))
\\
    \mu^{en}(\epsilon,c) &= \frac{\partial\psi}{\partial\epsilon}
    &= \mu_{ref} + \frac{R\theta_{ref}}{c_m}(c-c_{ref}) - \alpha^{ch} \colon \sigma(\epsilon,c)
\\
    j(\nabla[\mu])
    &:= -M \cdot \nabla[\mu]
\end{align}
```

## Weak Format

```math
\begin{align}
    \int_{\Omega} \sigma (\epsilon[u],c) : \epsilon[\delta u] ) \ d\Omega  &=  \int_{\Gamma_N^{(u)}} t^{p} \cdot \delta u \ d\Gamma   
    &
    \forall \delta u \in \mathbf{U}^{0}
\\
    \int_{\Omega} \dot{c} \delta \mu \ d\Omega - \int_{\Omega} j(\nabla[\mu]) \cdot \nabla[\delta \mu] \ d\Omega 
    &=  \int_{\Gamma_N^{(\mu)}} h^{p} \delta \mu d\Gamma 
    &
    \forall \delta \mu \in \mathbf{M}^{0}
\\
    \int_{\Omega} (\mu - \mu^{en}(\epsilon[u], c)) \delta c \ d\Omega 
    &= 0
    &
    \forall \delta \mu \in \mathbf{C}^{0}
\end{align}
```