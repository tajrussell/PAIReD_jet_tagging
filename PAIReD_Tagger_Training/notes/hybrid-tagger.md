# Details on the hybrid tagger implementation

The hybrid tagger implemented in the file [`PAIReD_ParT_sv_hybrid.py`](../networks/PAIReD_ParT_sv_hybrid.py) has the following details. 

#### Network and inputs
The architecture is given by the ParticleTransformer with both particle and SV input information.

#### Output nodes
The number of output nodes is given by:

$$N_\text{out} = N_\text{classes}+3\times N_\text{regression}$$

Hereby, $N_\text{classes}$ denotes the number of labels/classes for classification and $N_\text{regression}$ the number of variables to be regressed (e.g. the jet mass). The number of regressions is multiplied by three since there is one node for the regressed value $x_\text{reg}$ and two nodes for the error estimate (namely the $1\text{-}\sigma$ interval $[x_\text{err}^-, x_\text{err}^+]$) respectively.

#### Loss function
The implemented loss function $L$ is the sum of three loss terms:
```math
\boxed{L = L_\text{cls} + \lambda_\text{reg} \cdot L_\text{reg} + \lambda_\text{err} \cdot L_\text{err}}
```

Hereby, $L_\text{cls}$ is the loss for the classification, $L_\text{reg}$ the loss for the mass regression and $L_\text{err}$ the loss for the error estimation. $\lambda_\text{reg/err}$ are scaling factors that can be adjusted to set the balance between the different terms.

The cross-entropy loss is selected as the classification loss $L_\text{cls}$. This includes all jets equally.

In contrast, the regression-related losses make a strict distinction between resonant Higgs jets and non-resonant jets. This is because while there is a clear target value $x_\text{true}$ for the regression for Higgs jets (namely the MC-generated Higgs properties, e.g. $m_\text{H}$), non-resonant jets do not have a well-defined true value. Accordingly, the loss functions for non-resonant jets are set to baseline values $c_\text{reg/err}$ that have to be chosen. The regression loss $L_\text{reg}$ is then implemented with the log-cosh loss as follows:
```math
L_\text{reg} \;=\; \left\{ \begin{array}{ll}
    \log\left(\cosh\left(x_\text{true} - x_\text{reg} \right)\right), & \text{for Higgs PAIReD jets (CC or BB)} \\
    c_\text{reg}, & \text{otherwise}
\end{array} \right.
```

The loss for the error estimate is based on the quantile loss function $\rho_\tau$:
```math
\rho_\tau (\delta) \;=\; \left\{ \begin{array}{ll} \tau \delta, & \text{if }\delta >0 \\ (\tau -1)\,\delta, & \text{otherwise} \end{array} \right.
```

The total error loss $L_\text{err}$ is then given by the sum of the losses for the 16% and 84% quantiles:
```math
L_\text{err} \;=\; \left\{ \begin{array}{ll} \rho_{0.16}\left(x_\text{true} - x_\text{err}^-\right) \; + \; \rho_{0.84}\left(x_\text{true} - x_\text{err}^+\right), & \text{for Higgs PAIReD jets (CC or BB)} \\  c_\text{err}, & \text{otherwise} \end{array} \right.
```

Ideal values for $\lambda_\text{reg/err}$ and $c_\text{reg/err}$ will have to be determined when training on variable-mass samples by checking the converging values for the different loss terms. They can be set for training like this:
- $\lambda_\text{reg}$ is set by `--loss-option "factor_reg" $VALUE`
- $\lambda_\text{err}$ is set by `--loss-option "factor_err" $VALUE`
- $c_\text{reg}$ is set by `--loss-option "baseline_reg" $VALUE`
- $c_\text{err}$ is set by `--loss-option "baseline_err" $VALUE`