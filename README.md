WRI with Uncertainties
=========

This code provides the basic building blocks to test optimization algorithms on seismic inverse problems.
The canonical seismic waveform inversion problem is given by

$$\min_{m} \sum_{i} \frac{1}{2}||P^TA(m)^{-1}q_i - d_i||_2^2 + \frac{\alpha}{2}||Lm||_2^2,$$

where $A(m)$ is the discretized Helmholtz operator $\omega^2 m + \nabla^2$ with absorbing boundary conditions, 
$P$ is the sampling operator, $q_i$ are the source vectors, $d_i$ are the observed data and $L$ is the discretized 
$\nabla$.

The main function is *misfit.m*, which takes as input the medium parameters $m$ and the data $d_i$ (and definitions of $P$, $q_i$, $\omega$ etc.) and returns the misfit value, its gradient and a function handle to compute the action of the Gauss-Newton hessian.

For an example, see */examples/marm.m*, which uses a simple BB iteration to solve the optimization problem. 
Replace this with your favourite optimization algorithm and you're ready to roll!

Feel free to contact me with any questions, suggestions, etc.



Tristan van Leeuwen - T.vanLeeuwen@uu.nl





This code was used to produce the examples in the paper: Mitigating local minima in full-waveform inversion by expanding the search space, T. van Leeuwen and F.J. Herrmann, 2013 (submitted to GJI).

The code is distributed under the GNU GPL so that others may reproduce these results and conduct additional numerical experiments or modify the code to suit their needs.

To reproduce the examples from the paper, start Matlab and run example*.m.