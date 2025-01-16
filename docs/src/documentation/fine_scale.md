# [Chemo-Mechanical Problem - fine scale](@id documentation-fine-scale)

## Strong From

Three primary variable fields: displacement vector $\boldsymbol u(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}^{3}$
,molar ion concentration $\boldsymbol c(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$
as well as the chemical potential gradient $\boldsymbol \mu(x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$
are solved in the lithium structural battery problem.
To define the problem the balance of linear momentum and the balance of mass are given in Strong Forms as following:

```math
\begin{align}
-\boldsymbol\sigma \cdot \boldsymbol\nabla &= 0
&\
in \ \Omega \times (0,T]
\\
\dot{\boldsymbol c} + \boldsymbol j \cdot \boldsymbol\nabla &= 0
&\
in \ \Omega \times (0,T]
\end{align}
```

Where the boundary is divied into three parts for the three variable fields respectively
($\Gamma = \Gamma_D^{(u)} \cup \Gamma_N^{(u)} = \Gamma_D^{(\mu)} \cup \Gamma_N^{(\mu)} $):


```math
\begin{align}
\boldsymbol u &= \boldsymbol u^\text{p}
&\
on \ \Gamma_{D}^{(u)} \times (0,T]
\\
\boldsymbol t := \boldsymbol \sigma \cdot \boldsymbol n &= \boldsymbol t^\text{p}
&\
on \ \Gamma_{N}^{(u)} \times (0,T]
\\
\boldsymbol \mu &= \boldsymbol \mu^\text{p}
&\
on \ \Gamma_{D}^{(\mu)} \times (0,T]
\\
\boldsymbol h := -\boldsymbol j \cdot \boldsymbol n &= \boldsymbol h^\text{p}
&\
on \ \Gamma_{N}^{(\mu)} \times (0,T]
\end{align}
```

## Constitutive Models
Due to the coupling of the mechanical and chemical aspects the total strain has both the contribution due to a deformation state $\boldsymbol u(x,t)$
and the ion concentration $\boldsymbol c(x,t)$, where an ion intercalation tensor $\boldsymbol \alpha^{ch} = \boldsymbol \alpha \boldsymbol I$ is used based on the models of Zhang, Shyy and Sastry and Bohn et al.

```math
\begin{align}
\boldsymbol \varepsilon[\boldsymbol u] &:= (\boldsymbol u \otimes \boldsymbol \nabla)^\text{sym}
\\
\boldsymbol \varepsilon^\text{ch}(c) &:= \boldsymbol \alpha^\text{ch} (\boldsymbol c-\boldsymbol c_\text{ref})
\end{align}
```
The free energy $\psi$ is assumed to be the sum of mechanical and chemical parts
$\boldsymbol \psi (\boldsymbol \varepsilon, \boldsymbol c) = \boldsymbol \psi^\text{mech}(\boldsymbol \varepsilon,\boldsymbol c) + \boldsymbol \psi^\text{chem}(\boldsymbol c)$.

```math
\begin{align}
\boldsymbol \psi^\text{mech} &:= \frac{1}{2}(\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}(\boldsymbol c)) \colon \boldsymbol E \colon (\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}(\boldsymbol c))
\\
\boldsymbol \psi^\text{chem} &:= (\boldsymbol c-\boldsymbol c_\text{ref})\boldsymbol \mu_\text{ref} + \frac{1}{2} \frac{\boldsymbol R \boldsymbol \theta_\text{ref}}{\boldsymbol c_\text m}(\boldsymbol c-\boldsymbol c_\text{ref})^2
\end{align}
```

Furthermore, a constant mobility tensor $\boldsymbol M$ for the assumption of a linear relation between ion flux and the gradient of
the chemical potential is defined using a mobility coefficient $\boldsymbol \eta$ for the isotropic case such as in this project.
A reference temperature $\boldsymbol \theta_\text{ref}$ and the concentration $\boldsymbol c_\text{ref}$ as well as the reference chemical potential $\boldsymbol \mu_\text{ref}$ are gloable constant material parameters for the purpose of linearization. 
So that the simplified constitutive equations are:

```math
\begin{align}
    \boldsymbol \sigma(\boldsymbol \varepsilon,\boldsymbol c) &= \frac{\partial \boldsymbol\psi}{\partial \boldsymbol\varepsilon} = \boldsymbol E \colon (\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}(\boldsymbol c))
\\
    \boldsymbol \mu(\boldsymbol \varepsilon,\boldsymbol c) &= \boldsymbol \mu^\text{en}(\boldsymbol \varepsilon,\boldsymbol c) = \frac{\partial \boldsymbol \psi}{\partial \boldsymbol \varepsilon} = \boldsymbol \mu_\text{ref} + \frac{\boldsymbol R \boldsymbol \theta_\text{ref}}{\boldsymbol c_m}(\boldsymbol c-\boldsymbol c_\text{ref}) - \boldsymbol \alpha^\text{ch} \colon \boldsymbol \sigma(\boldsymbol \varepsilon,\boldsymbol c)
\\
    \boldsymbol j(\boldsymbol \nabla[\boldsymbol \mu]) &:= -\boldsymbol M \cdot \boldsymbol \nabla[\boldsymbol \mu]
\end{align}
```

## Weak Format

```math
\begin{align}
\int_{\Omega} \boldsymbol \sigma (\boldsymbol \varepsilon[\boldsymbol u],\boldsymbol c) : \boldsymbol \varepsilon[\delta \boldsymbol u] \ d\Omega  &=  \int_{\Gamma_N^{(u)}} \boldsymbol t^\text{p} \cdot \delta \boldsymbol u \ d\Gamma
&\
\forall \delta \boldsymbol u \in \mathbb{U}^{0}
\\
\int_{\Omega} \dot{\boldsymbol c} \ \delta \boldsymbol \mu \ d\Omega - \int_{\Omega} \boldsymbol j(\boldsymbol \nabla[\boldsymbol \mu]) \cdot \boldsymbol \nabla[\delta \boldsymbol \mu] \ d\Omega
&=  \int_{\Gamma_N^{(\mu)}} \boldsymbol h^\text{p} \delta \boldsymbol \mu \ d\Gamma
&\
\forall \delta \boldsymbol \mu \in \mathbb{M}^{0}
\\
\int_{\Omega} (\boldsymbol \mu - \boldsymbol \mu^\text{en}(\boldsymbol \varepsilon[u], \boldsymbol c)) \delta \boldsymbol c \ d\Omega
&= 0
&\
\forall \delta \boldsymbol \mu \in \mathbb{C}^{0}
\end{align}
```
where the initial state of all unknown fields for all elements in $\Omega$ are choosen based on the consistent reference state as mentioned earlier $\boldsymbol c_0 = \boldsymbol c_\text{ref}, \boldsymbol \mu_0 = \boldsymbol \mu_\text{ref}, \boldsymbol u_0 = 0$.

## Time stepping

As the transit problem is discussed here, one of the crucial part is to solve the time derivative term $\int_{\Omega} \dot{\boldsymbol c} \ \delta \boldsymbol \mu \ d\Omega$ in the FE context. Due to the linearity of the problem a direct derivation of the element stiffness $\boldsymbol K_e$, the mass matrix $\boldsymbol M_e$, and the right hand side vector  $\boldsymbol f_e$ is possible. The constraint matrix $\boldsymbol C_e$ for $\langle \boldsymbol \mu\rangle _{\square} = \bar{\boldsymbol \mu}$ is merged into the last row and colum of $\boldsymbol K_e$. 

The whole PDE system can be discretized using the Crank-Nicolson scheme, where $\boldsymbol a$ is the global result vector including all unknown fields:
```math
\begin{align}
\boldsymbol a(t^\text{n+1}) &\approx \boldsymbol a(t^\text{n}) + \frac{1}{2} \Delta t (\dot {\boldsymbol a} (t^{n}) + \dot {\boldsymbol a}(t^{n+1}))
\end{align}
```
so that the linear system can be discretized as following:
```math
\begin{align}
\boldsymbol M \dot{\boldsymbol a} + \boldsymbol K \boldsymbol a^\text{n} &= \boldsymbol f
\\
\dot{\boldsymbol a} &= \boldsymbol M^{-1}[\boldsymbol f - \boldsymbol K \boldsymbol a^\text{n}]
\\
\boldsymbol a^\text {n+1} &= \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (\dot {\boldsymbol a} ^\text {n} + \dot {\boldsymbol a}^\text {n+1})
\\
\boldsymbol a^\text {n+1} &= \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (\boldsymbol M^{-1}(2 \boldsymbol f  - \boldsymbol K (\boldsymbol a^\text{n} + \boldsymbol a^\text{n+1})) )
\\
\boldsymbol M \boldsymbol a^\text {n+1} + \frac{1}{2} \Delta t \boldsymbol K \boldsymbol a^\text{n+1} &= \boldsymbol M \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (-\boldsymbol K \boldsymbol a^\text{n} + 2\boldsymbol f)
\end{align}
```
In order to solve for $\boldsymbol a^\text{n+1}$ the jacobian matrix $\boldsymbol J$ and the residual vector $\boldsymbol g$ are computed using the gloable stiffness matrix $\boldsymbol K$ and mass matrix $\boldsymbol M$ as well as $\boldsymbol a^\text{n}$ for $\boldsymbol g$.
```math
\begin{align}
\boldsymbol J &= \boldsymbol M + \frac{1}{2}\Delta t \boldsymbol K
\\
\boldsymbol g &= (\boldsymbol M -\frac{1}{2}\Delta t \boldsymbol K )\boldsymbol a^\text{n} + \Delta t \boldsymbol f
\end{align}
```
where the initial state of all unknown fields for all elements in $\Omega$ are choosen as $c_0 = c_{ref}, \mu_0 = \mu_{ref}, u_0 = 0$.
