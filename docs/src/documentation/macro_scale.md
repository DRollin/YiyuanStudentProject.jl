# [Chemo-Mechanical Problem - macro scale](@id documentation-macro-scale)
## Weak Format
```math
\begin{align}
\int_{\Omega} \bar{\boldsymbol \sigma} : \boldsymbol \varepsilon[\delta \bar{\boldsymbol u}] \ d\Omega  &=  \int_{\Gamma_N^{(u)}} \boldsymbol t^\text{p} \cdot \delta \bar{\boldsymbol u} \ d\Gamma
&\
\forall \delta \bar{\boldsymbol u} \in \bar{\mathbb{U}}^{0}
\\
\int_{\Omega} \dot{\bar{\boldsymbol c}} \ \delta \bar{\boldsymbol \mu} \ d\Omega - \int_{\Omega} (\bar{\boldsymbol j} - \dot{\bar{\boldsymbol c}}_{2}) \cdot \boldsymbol \nabla[\delta \bar{\boldsymbol \mu}] \ d\Omega
&=  \int_{\Gamma_N^{(\mu)}} \boldsymbol h^\text{p} \delta \bar{\boldsymbol \mu} \ d\Gamma
&\
\forall \delta \bar{\boldsymbol \mu} \in \bar{\mathbb{M}}^{0}
\end{align}
```

with the variationally consistent macro scale homogenized fields:

```math
\begin{align}
    \bar{\boldsymbol \sigma} &:= \lang \boldsymbol \sigma \rang _\square = \lang \boldsymbol E \colon [\boldsymbol \varepsilon[\boldsymbol u]-\boldsymbol \alpha^\text{ch}[\boldsymbol c - c_{ref}]] \rang _\square
\\
    \bar{\boldsymbol j} &:= \lang \boldsymbol j \rang _\square = -\lang \boldsymbol M \cdot \boldsymbol \zeta[\boldsymbol \mu] \rang _\square
\\
    \bar{\boldsymbol c} &:= \lang \boldsymbol c \rang _\square
\\
    \bar{\boldsymbol c_{2}} &:= \lang \boldsymbol c [\boldsymbol x - \bar{\boldsymbol x}] \rang _\square
\end{align}
```
## Time Stepping
For every time step the RVE problem is going to be solved at the corresponding quadrature point. The four variationally consistent macro scale fields are updated using new RVE results at current time step. New stiffness matrix $\boldsymbol K^{n}$ is then computed with the new variationally consistent macro scale fields. A time dependent boundary condition is applied on both $\boldsymbol K^{n}$ and $\boldsymbol f^{n}$.

```math
\begin{align}
\boldsymbol K^{n} \boldsymbol a^\text{n} &= \boldsymbol f^{n}
\\
\boldsymbol a^{n} &= (\boldsymbol K^{n})^{-1} \boldsymbol f^{n}

\end{align}
```