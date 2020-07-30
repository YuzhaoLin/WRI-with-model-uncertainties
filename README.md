WRI with Uncertainties
=========

This code provides the basic building blocks to test optimization algorithms on seismic inverse problems.
The canonical seismic waveform inversion problem is given by

$$ \phi (\mathbf{m}) =
    \left \|
         P_i A (\mathbf{m})^{-1} \mathbf{q}_i - \mathbf{d}_i
    \right \| ^{2} _{\Sigma(\mathbf{m})_i + \Sigma_{m}}$$

where $A(m)$ is the discretized Helmholtz operator $\omega^2 m + \nabla^2$ with absorbing boundary conditions, 
$P$ is the sampling operator, $q_i$ are the source vectors, $d_i$ are the observed data and $L$ is the discretized 
$\nabla$. And $\Sigma(m)$ is a weight function consists by theory covariance matrix:

$$\Sigma(\mathbf{m})_i  = P_i( A( \textbf{m})^{*}
        \Sigma^{-1}_{p} A( \mathbf{m}))^{-1} P_i^{*}$$



The main function is *misfit.m*, which takes as input the medium parameters $m$ and the data $d_i$ (and definitions of $P$, $q_i$, $\omega$ etc.) and returns the misfit value, and its gradient. Three files corresponding three $\Sigma(m)$ forms.



This code was used to produce the basic examples in the paper: Wavefield reconstruction inversion with theory uncertainties, Yuzhao Lin and T. van Leeuwen, 2020 (submitted to GJI).

The code is distributed under the GNU GPL so that others may reproduce these results and conduct additional numerical experiments or modify the code to suit their needs.



Feel free to contact me with any questions, suggestions, etc.

Yuzhao Lin - y.lin2@uu.nl

Tristan van Leeuwen - T.vanLeeuwen@uu.nl





