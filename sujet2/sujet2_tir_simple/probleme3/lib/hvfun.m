function hv = hvfun(t, z, par)
%-------------------------------------------------------------------------------------------
%
%    hvfun
%
%    Description
%
%        Computes the Hamiltonian vector field associated to H.
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        hv = hvfun(t, z, par)
%
%    Inputs
%
%        t    -  real        : time
%        z    -  real vector : state and costate
%        par  -  real vector : parameters, par=[] if no parameters
%
%    Outputs
%
%        hv   -  real vector : hamiltonian vector field at time t
%
%-------------------------------------------------------------------------------------------
n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);
%% A REMPLACER
u = control(t, z, par);
hv = [x((length(x)/2)+1:length(x)); u; 0.0 ; -p(1:(length(p)/2))];

