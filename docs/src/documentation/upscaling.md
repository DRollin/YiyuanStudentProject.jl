# [Chemo-Mechanical Problem - upscaling](@id documentation-upscaling)
## Running Averaging on RVE
The homogenization of macro scale problem is introduced by using the running averages on RVE. This approach for a $\diamond$ quantity is denoted for volume and surface as following respectively:

```math
\begin{align}
\lang \diamond \rang _{\square} &:= \frac{1}{|\boldsymbol \Omega _{\square}|} \int_{{\boldsymbol \Omega} _{\square}} \diamond \ d \boldsymbol \Omega

\\
\lang \lang \diamond \rang \rang _{\square} &:= \frac{1}{|\boldsymbol \Gamma _{\square}|} \int_{{\boldsymbol \Gamma} _{\square}} \diamond \ d \boldsymbol \Gamma
\end{align}
```
where the corresponding homogenized macro scale field $\bar \diamond$ is thus defined as $\bar \diamond := \lang \diamond \rang _{\square}$ . This approach allows the fine(sub) scale quantities upscale into the macro scale.

## Hierarchical Decomposition
For the chemo-mechanical problem, unknown fields displacement $\boldsymbol u$ and chemical potential $\boldsymbol \mu$ can be decomposed into (smooth) macro-scale $\diamond ^{M}$ and (fluctuating) sub-scale $\diamond ^{s}$, where the macro-scale parts will be prescribed using Taylor expensions in each RVE. Since the sub-scale unknown field concentration $\boldsymbol c$ will not be considered in macro scale problem due to its localization, no such decomposition is needed for it.

Thus, using first order Taylor expension(linear variation) the homogenized fields $\bar{\boldsymbol u}$ and $\bar{\boldsymbol \mu}$ are defined as follows:

```math
\begin{align}
\boldsymbol u^{M}[\bar{\boldsymbol u}](\bar{\boldsymbol x}, \boldsymbol x ) &= \bar{\boldsymbol u}(\bar{\boldsymbol x}) + \bar{\boldsymbol \xi}(\bar{\boldsymbol x}) \cdot [x - \bar{\boldsymbol x}], \ \ \bar{\boldsymbol \xi} := \bar{\boldsymbol u}  \otimes \boldsymbol \nabla, & \ x \in \boldsymbol \Omega_\square

\\

\boldsymbol \mu^{M}[\bar{\boldsymbol \mu}](\bar{\boldsymbol x}, \boldsymbol x ) &= \bar{\boldsymbol \mu}(\bar{\boldsymbol x}) + \bar{\boldsymbol \zeta}(\bar{\boldsymbol x}) \cdot [x - \bar{\boldsymbol x}], \ \ \bar{\boldsymbol \zeta} := \boldsymbol \nabla \bar{\boldsymbol \mu}  , & \ x \in \boldsymbol \Omega_\square
\end{align}
```
where $\bar{\boldsymbol x}$ is centriod and also holds to the upscaling rules mentioned earlier.

## Upscaling
After inserting the hierarchical decomposition $\boldsymbol u = \boldsymbol u ^{M} + \boldsymbol u ^{s} \in \mathbb U ^{M} \oplus \mathbb U ^{s}$ and $\boldsymbol \mu = \boldsymbol \mu ^{M} + \boldsymbol \mu ^{s} \in \mathbb M ^{M} \oplus \mathbb M ^{s}$ and the upsacling assumption $\bar \diamond := \lang \diamond \rang _{\square}$ into the fine scale week form:

```math
\begin{align}
\int_{\Omega} \lang \boldsymbol \sigma (\boldsymbol \varepsilon[\boldsymbol u ^{M}[\bar {\boldsymbol u}] + \boldsymbol u ^{s}],\boldsymbol c ^{s}) : \boldsymbol \varepsilon[\boldsymbol u ^{M}[ \delta \bar {\boldsymbol u} ]+ \delta\boldsymbol u ^{s}] \rang _\square \ d\Omega  &=  \int_{\Gamma_N^{(u)}} \boldsymbol t^\text{p} \cdot \delta \bar{\boldsymbol u} \ d\Gamma
&\
\forall \delta (\bar{\boldsymbol u}, {\boldsymbol u} ^{s}) &\in \bar {\mathbb{U}}^{0} \oplus \mathbb U ^{s}
\\
\int_{\Omega} \lang \dot{\boldsymbol c} ^{s} \ [\mu ^{M} [\delta \bar{\boldsymbol \mu}] + \delta\mu^{s}] \rang _\square \ d\Omega - \int_{\Omega} \boldsymbol \lang j(\boldsymbol \nabla[\mu ^{M} [\bar{\boldsymbol \mu}] + \mu^{s}]) \cdot \boldsymbol \nabla[\mu ^{M} [\delta \bar{\boldsymbol \mu}] + \delta\mu^{s}] \rang _\square \ d\Omega
&=  \int_{\Gamma_N^{(\mu)}} \boldsymbol h^\text{p} \delta \bar{\boldsymbol \mu} \ d\Gamma
&\
\forall \delta (\bar{\boldsymbol \mu}, {\boldsymbol \mu} ^{s}) &\in \bar {\mathbb{M}}^{0} \oplus \mathbb M ^{s}
\\
\int_{\Omega} \lang ([\mu ^{M} [\bar{\boldsymbol \mu}] + \mu^{s}] - \boldsymbol \mu^\text{en}(\boldsymbol \varepsilon[\boldsymbol u ^{M}[\bar {\boldsymbol u}] + \boldsymbol u ^{s}], \boldsymbol c^{s})) \delta \boldsymbol c^{s} \rang _\square \ d\Omega
&= 0
&\
\forall \delta \boldsymbol c^{s} &\in \mathbb{C}^{0}
\end{align}
```
so that with setting $\delta \boldsymbol u^{s}, \delta \boldsymbol \mu^{s}, \delta \boldsymbol \mu^{s}$ to 0, the globally coupled macro scale problem is defined.