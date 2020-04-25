function du = dcontrol(t, z, par)
%-------------------------------------------------------------------------------------------
%
%    dcontrol
%
%    Description
%
%        Computes the derivative of the control.
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        du = dcontrol(t, z, par)
%
%    Inputs
%
%        t    -  real        : time
%        z    -  real vector : state and costate
%        par  -  real vector : parameters, par=[] if no parameters
%
%    Outputs
%
%        du   -  real vector : the derivative of the control
%
%-------------------------------------------------------------------------------------------

% par = [t0; tf; x0; xf, epsilon]
%
n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);

%% A REMPLACER
du = 0.0;
