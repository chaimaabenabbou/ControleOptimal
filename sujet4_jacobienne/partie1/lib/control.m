function u = control(t, z, par)
%-------------------------------------------------------------------------------------------
%
%    control
%
%    Description
%
%        Computes the control.
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        u = control(t, z, par)
%
%    Inputs
%
%        t    -  real        : time
%        z    -  real vector : state and costate
%        par  -  real vector : parameters, par=[] if no parameters
%
%    Outputs
%
%        u   -  real vector : control
%
%-------------------------------------------------------------------------------------------

% par = [t0; tf; x0; xf, epsilon]
%

n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);
eps0 = par(5) ;
psip = -1+abs(p);
if p ~= zeros(size(p))
    u = -2*eps0*sign(p)/(psip-2*eps0-sqrt(psip^2+4*eps0^2));
else
    u = -2*eps0/(-1-2*eps0-sqrt(1+4*eps0^2));
end