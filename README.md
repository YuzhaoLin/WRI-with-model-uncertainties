#####  

In seismic forward modeling theory, under the assumption of the acoustic medium, wave propagation can be described by the acoustic wave equation:
$$
\mathbf{A(\mathbf{m})}
\mathbf{u}=
\mathbf{q}
$$
in which $\mathbf{u}$ is wavefield, $\mathbf{A}=\omega^2diag(\mathbf{m})+L$ represent a frequency-domain forward modeling operator, $\mathbf{m}$ is model parameter, $\mathbf{q}$ denotes source and $L$ refers to the discretized Laplacian.

Weglein (2013) pointed out that in the application of full waveform inversion, there are some errors in the data and the earth model. To quantitatively describe such problems, the observed data can be written as:
$$
\mathbf{d}=
\mathbf{P}
(\mathbf{A}(\mathbf{m}) + \epsilon_p )^{-1}
\mathbf{q} 
+\epsilon_m
$$
in which $\mathbf{d}$ is the observed data and $\mathbf{P}$ is the sampling operator.  $\epsilon_p$ , $\epsilon_m$ are used to describe the theoretical uncertainties and the measurement uncertainties , which satisfies:
$$
\epsilon_m \sim
\mathscr{N}(0,\Sigma_m)
$$

$$
\epsilon_p \sim 
\mathscr{N}(0,\Sigma_p)
$$

where $\mathscr{N}$ represent the Gaussian distribution and $\Sigma$ is the corresponding covariance matrix. 



The traditional FWI method (Tarantola, 1984; Pratt, 1999) is to minimize a nonlinear least-squares problem, which ignores theoretical and measurement uncertainties:
$$
\underset{ \mathbf{m} }  { min }
\left \| 
\mathbf{P}\mathbf{A}( \mathbf{m})^{-1} \mathbf{q}- \mathbf{d}
\right \| ^{2}
$$

Obviously, equation (5) is unable to provide reliable inversion results when the covariance matrices seriously affects forward modeling and measurement. 



##### WRI method

In order to accurately describe the inversion problem with uncertainties, a wavefield reconstruction inversion (WRI) objective function can be obtained (van Leeuwen and Herrmann, 2013, 2015) based on the Bayesian interference equation and Gaussian distribution assumption (Tarantola, 2005) (see Appendix A for details):
$$
\underset{ \mathbf{u}, \mathbf{m} }  { min }
\left \| 
\mathbf{P}\mathbf{u} - \mathbf{d}
\right \| ^{2} _{\Sigma_{m}} +
\left \| 
\mathbf{A}( \mathbf{m}) \mathbf{u} - \mathbf{q}
\right \|^{2}_{\Sigma_{p}}
$$
##### 

In this section, an alternative (but equivalent) formula for extended full-waveform inversion from WRI with uncertainties is derived. This reformulation takes the form of a conventional FWI formula with a medium-dependent weight, which is more easy and cheap to compute.

Firstly, based on the WRI method, considering measurement and physical uncertainties, the extended source formulation can be obtained by introducing a new source-like variable:
$$
\mathbf{r} = 
\mathbf{A}( \mathbf{m}) \mathbf{u} - \mathbf{q}
$$
The inversion problem of equation (6) becomes:
$$
\underset{ \mathbf{r}, \mathbf{m} }  { min }
\left \| 
\mathbf{P}\mathbf{A}(\mathbf{m})^{-1} (\mathbf{r} + \mathbf{q}) - \mathbf{d}
\right \| ^{2} _{\Sigma_{m}} +
\left \| 
\mathbf{r}
\right \|^{2}_{\Sigma_{p}}
$$

Minimizing the new variable gives:
$$
((\mathbf{P}\mathbf{A}(\mathbf{m})^{-1})^{*}
\Sigma_{m}^{-1}
(\mathbf{P}\mathbf{A}(\mathbf{m})^{-1})
+\Sigma_{p}^{-1}  )
\mathbf{r}_0 
= ((\mathbf{P}\mathbf{A}(\mathbf{m})^{-1})^{*}
\Sigma_{m}^{-1}
(\mathbf{d}-\mathbf{P}\mathbf{A}(\mathbf{m})^{-1}\mathbf{q}))
$$

which is the traditional form of the extended FWI. Note that different $\Sigma_p$ choices lead to different extended FWI methods. For example, let $\Sigma_p^{-1}$ be a function of source location distance in the frequency-domain, we will have the source extended FWI (Huang et al., 2018). However, the covariance matrices in equation (11) are independent variables and need to be defined according to our theoretical knowledge.



Via some algebra(see Appendix B for more details), we can obtain an objection function in the traditional FWI form with a medium-dependent weight:
$$
\phi (\mathbf{m}) = 
\left \| 
\mathbf{P}\mathbf{A}(\mathbf{m})^{-1} \mathbf{q} - \mathbf{d}
\right \| ^{2} _{\Sigma(\mathbf{m})+\Sigma_{m}}
$$
with:
$$
\Sigma(\mathbf{m}) 
=
\mathbf{P}(
\mathbf{A}( \textbf{m})^{*} 
\Sigma^{-1}_{p} 
\mathbf{A}( \mathbf{m}))^{-1}
\mathbf{P}^{*}
$$

In the new objective function, one can see that the theoretical covariance matrix and the measurement covariance matrix constitute a weight function. 



Remarkably, the gradient of the objective has a simple expression, similar to that of conventional FWI:
$$
g(\mathbf{m})=
\frac { \partial \phi (\mathbf{m}) } { \partial \mathbf{m} } = 
\frac { \partial } { \partial \mathbf{m} }
( \mathbf{P}\mathbf{A}(\mathbf{m})^{-1} \mathbf{q} - \mathbf{d} )^*
( \Sigma(\mathbf{m})+\Sigma_{m} )^{-1}
( \mathbf{P}\mathbf{A}(\mathbf{m})^{-1} \mathbf{q} - \mathbf{d} )
$$
more specifically:
$$
g(\mathbf{m})=
-2\mathbf{u}^{*}_{0}  
\frac {\partial \mathbf{A}} {\partial \mathbf{m}}
\mathbf{v}_0  +
2\mathbf{v}^{*}_{0}  
\frac {\partial \mathbf{A}} {\partial \mathbf{m}}
\mathbf{w}_0
$$
in which:
$$
\mathbf{u}_{0}  = \mathbf{A}^{-1} \mathbf{q} 
$$

$$
\mathbf{w}_{0}  = \mathbf{A}^{-1} \Sigma_{p}  \mathbf{v}_0
$$

$$
\mathbf{v}_{0}  = ( \mathbf{P} \mathbf{A}^{-1})^* 
(\Sigma(\mathbf{m}) + \Sigma_{m}) ^{-1}
\mathbf{h}  = 
( \mathbf{P} \mathbf{A}^{-1})^* 
\widehat{\mathbf{h}}
$$

Note that the gradient can be readily computed using the standard tools available in any FWI workflow which only requires one additional forward modeling. The specific algorithm of the inversion will be given after the estimation of covariance matrices and the new residual. 


