# [Chemo-Mechanical Problem - Fine-Scale](@id documentation-fine-scale)

## Strong From

Three primary variable fields: displacement vector $\boldsymbol u(\boldsymbol x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}^{3}$, 
molar ion concentration $c(\boldsymbol x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$
as well as the chemical potential $ \mu(\boldsymbol x,t): \Omega \times \mathbb{R}^{+} \rightarrow \mathbb{R}$
are considered.
To define the problem the balance of linear momentum and the balance of mass are given in strong forms as:

```math
\begin{align}
-\boldsymbol\sigma \cdot \boldsymbol\nabla &= 0
&\
\text{in} \ \Omega \times (0,T]
\\
\dot{ c} + \boldsymbol j \cdot \boldsymbol\nabla &= 0
&\
\text{in} \ \Omega \times (0,T]
\end{align}
```

Where the boundary is divied into three parts for the three variable fields respectively
$\Gamma = \Gamma_D^{(u)} \cup \Gamma_N^{(u)} = \Gamma_D^{(\mu)} \cup \Gamma_N^{(\mu)}$:


```math
\begin{align}
\boldsymbol u &= \boldsymbol u^\text{p}
&\
\text{on} \ \Gamma_{D}^{(u)} \times (0,T]
\\
\boldsymbol t := \boldsymbol \sigma \cdot \boldsymbol n &= \boldsymbol t^\text{p}
&\
\text{on} \ \Gamma_{N}^{(u)} \times (0,T]
\\
 \mu &=  \mu^\text{p}
&\
\text{on} \ \Gamma_{D}^{(\mu)} \times (0,T]
\\
\boldsymbol h := -\boldsymbol j \cdot \boldsymbol n &= \boldsymbol h^\text{p}
&\
\text{on} \ \Gamma_{N}^{(\mu)} \times (0,T]
\end{align}
```

## Constitutive Models
Due to the coupling of the mechanical and chemical processes the total strain contains contributions due to a deformation state $\boldsymbol u(\boldsymbol x,t)$
and the ion concentration $c(\boldsymbol x,t)$, where an ion intercalation tensor $\boldsymbol \alpha^{ch} = \alpha \boldsymbol I$ is used based on the models of [Zhang.2007](refs.bib), [Bohn.2013](refs.bib).

```math
\begin{align}
\boldsymbol \varepsilon[\boldsymbol u] &:= (\boldsymbol u \otimes \boldsymbol \nabla)^\text{sym}
\\
\boldsymbol \varepsilon^\text{ch}(c) &:= \boldsymbol \alpha^\text{ch} ( c- c_\text{ref})
\end{align}
```
The free energy $\psi$ is assumed to be the sum of mechanical and chemical parts
$\psi (\boldsymbol \varepsilon, c) =  \psi^\text{mech}(\boldsymbol \varepsilon, c) +  \psi^\text{chem}( c)$.

```math
\begin{align}
 \psi^\text{mech} &:= \frac{1}{2}(\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}( c)) \colon \boldsymbol E \colon (\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}( c))
\\
 \psi^\text{chem} &:= ( c- c_\text{ref}) \mu_\text{ref} + \frac{1}{2} \frac{ R  \theta_\text{ref}}{ c_\text m}( c- c_\text{ref})^2
\end{align}
```

Furthermore, a constant mobility tensor $\boldsymbol M$ for the assumption of a linear relation between ion flux and the gradient of
the chemical potential is defined. For isotropic cases as in this project a mobility coefficient $\eta$ is used.
A reference temperature $\theta_\text{ref}$ and the concentration $c_\text{ref}$ as well as the reference chemical potential $\mu_\text{ref}$ are gloable constant material parameters for the purpose of linearization. $c_\text{m}$ is an auxiliary parameter.
So that the simplified constitutive equations are:

```math
\begin{align}
    \boldsymbol \sigma(\boldsymbol \varepsilon, c) &= \frac{\partial \boldsymbol\psi}{\partial \boldsymbol\varepsilon} = \boldsymbol E \ \colon (\boldsymbol \varepsilon-\boldsymbol \varepsilon^\text{ch}(c))
\\
    \mu(\boldsymbol \varepsilon, c) &=  \mu^\text{en}(\boldsymbol \varepsilon, c) = \frac{\partial \boldsymbol \psi}{\partial \boldsymbol \varepsilon} =  \mu_\text{ref} + \frac{ R  \theta_\text{ref}}{ c_\text{m}}( c- c_\text{ref}) - \boldsymbol \alpha^\text{ch} \colon \boldsymbol \sigma(\boldsymbol \varepsilon, c)
\\
    \boldsymbol j(\boldsymbol \nabla[ \mu]) &:= -\boldsymbol M \cdot \boldsymbol \nabla[ \mu]
\end{align}
```

## Weak Format

```math
\begin{align}
\int_{\Omega} \boldsymbol \sigma (\boldsymbol \varepsilon[\boldsymbol u], c) : \boldsymbol \varepsilon[\delta \boldsymbol u] \ d\Omega  &=  \int_{\Gamma_N^{(u)}} \boldsymbol t^\text{p} \cdot \delta \boldsymbol u \ d\Gamma
&\
\forall \delta \boldsymbol u \in \mathbb{U}^{0}
\\
\int_{\Omega} \dot{ c} \ \delta  \mu \ d\Omega - \int_{\Omega} \boldsymbol j(\boldsymbol \nabla[ \mu]) \cdot \boldsymbol \nabla[\delta  \mu] \ d\Omega
&=  \int_{\Gamma_N^{(\mu)}} \boldsymbol h^\text{p} \delta  \mu \ d\Gamma
&\
\forall \delta  \mu \in \mathbb{M}^{0}
\\
\int_{\Omega} ( \mu -  \mu^\text{en}(\boldsymbol \varepsilon[u],  c)) \delta  c \ d\Omega
&= 0
&\
\forall \delta  c \in \mathbb{C}^{0}
\end{align}
```
where the initial state $\diamond_{0}$ of all unknown fields for all elements in $\Omega$ can be choosen based on the consistent reference state as mentioned earlier $c_0 =  c_\text{ref},  \mu_0 = \mu_\text{ref},  u_0 = 0$.

## Time Stepping

As the transient problem is discussed here, one of the crucial part is to introduce a time integation $\int_{\Omega} \dot{ c} \ \delta  \mu \ d\Omega$ in the FE context. Due to the linearity of the problem, a direct derivation of the element stiffness $\boldsymbol K_\text{e}$, the mass matrix $\boldsymbol M_\text{e}$, and the right hand side vector  $\boldsymbol f_\text{e}$ is possible. The constraint matrix $\boldsymbol C$ for $\langle \boldsymbol \mu\rangle _{\square} = \bar{ \mu}$ is merged into the last row and column of $\boldsymbol K$. 

The whole PDE system can be discretized using the Crank-Nicolson scheme, where $\boldsymbol a$ is the global solution vector including all degrees of freedom:
```math
\begin{align}
\boldsymbol a(t^\text{n+1}) &\approx \boldsymbol a(t^\text{n}) + \frac{1}{2} \Delta t (\dot {\boldsymbol a} (t^{n}) + \dot {\boldsymbol a}(t^{n+1}))
\end{align}
```
so that the linear system can be discretized as follows:
```math
\begin{align}
\boldsymbol M \dot{\boldsymbol a^\text{n}} + \boldsymbol K \boldsymbol a^\text{n} &= \boldsymbol f^\text{n}
\\
\dot{\boldsymbol a^\text{n}} &= \boldsymbol M^{-1}[\boldsymbol f^\text{n} - \boldsymbol K \boldsymbol a^\text{n}]
\\
\boldsymbol a^\text {n+1} &= \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (\dot {\boldsymbol a} ^\text {n} + \dot {\boldsymbol a}^\text {n+1})
\\
\boldsymbol a^\text {n+1} &= \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (\boldsymbol M^{-1}((\boldsymbol f^\text{n} + \boldsymbol f^\text{n+1})  - \boldsymbol K (\boldsymbol a^\text{n} + \boldsymbol a^\text{n+1})) )
\\
\boldsymbol M \boldsymbol a^\text {n+1} + \frac{1}{2} \Delta t \boldsymbol K \boldsymbol a^\text{n+1} &= \boldsymbol M \boldsymbol a^\text {n} + \frac{1}{2} \Delta t (-\boldsymbol K \boldsymbol a^\text{n} + (\boldsymbol f^\text{n} + \boldsymbol f^\text{n+1}))
\end{align}
```
In order to solve for $\boldsymbol a^\text{n+1}$ the Jacobian matrix $\boldsymbol J$ and the residual vector $\boldsymbol g$ are computed using the global stiffness matrix $\boldsymbol K$ and mass matrix $\boldsymbol M$ as well as $\boldsymbol a^\text{n}$ for $\boldsymbol g$.
```math
\begin{align}
\boldsymbol J &= \boldsymbol M + \frac{1}{2}\Delta t \boldsymbol K
\\
\boldsymbol g &= (\boldsymbol M -\frac{1}{2}\Delta t \boldsymbol K )\boldsymbol a^\text{n} + \Delta t \frac{1}{2}(\boldsymbol f^\text{n} + \boldsymbol f^\text{n+1})
\\
\boldsymbol J \boldsymbol a^\text {n+1} &= \boldsymbol g
\end{align}
```
where the global right hand side vector $\boldsymbol f$ is time independent due to the assumption of fixed macro scale data on RVE. It only consists of constant material parameters like $c _\text{ref}, R, \theta_\text{ref}, c_m, \mu _\text{ref}, \boldsymbol E, \boldsymbol \alpha^\text{ch}$ in each element right hand side vector $\boldsymbol f_\text{e}$

```math
\begin{align}
{f_\text{e}}_{u} &= \int_{\Omega} - \delta \boldsymbol \varepsilon : \boldsymbol E : \boldsymbol \alpha^\text{ch} c _\text{ref} \ d\Omega  
\\
{f_\text{e}}_{c} &= \int_{\Omega} \delta c (\mu _\text{ref} - (\frac{ R  \theta_\text{ref}}{ c_\text{m}} + \boldsymbol \alpha^\text{ch} : \boldsymbol E : \boldsymbol \alpha^\text{ch} ))   c _\text{ref} \ d\Omega
\end{align}
```
