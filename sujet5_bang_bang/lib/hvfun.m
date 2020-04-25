function hv = hvfun(t, z, par, iarc)
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
%        hv = hvfun(t, z, par, iarc)
%
%    Inputs
%
%        t    -  real        : time
%        z    -  real vector : state and costate
%        par  -  real vector : parameters, par=[] if no parameters
%        iarc -  integer     : index of current arc
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
 u = control(t, z, par,iarc);
hv = [-x+u ; p];
