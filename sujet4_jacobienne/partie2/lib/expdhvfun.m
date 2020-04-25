function [ tout, z, dz ]  = expdhvfun(tspan, z0, dz0, options, par, derivativeChoice)
%-------------------------------------------------------------------------------------------
%
%    exphvfun
%
%    Description
%
%        Computes the chronological exponential of the variational system associated to hv.
%
%-------------------------------------------------------------------------------------------
%
%    Matlab Usage
%
%        [tout, z, dz, flag] = expdhvfun(tspan, z0, dz0, options, par, choix)
%
%    Inputs
%
%        tspan   - real row vector of dimension 2   : tspan = [t0 tf]
%        z0      - real vector                      : initial flow
%        dz0     - real vector                      : initial value of the variational eq
%        options - struct vector                    : odeset options
%        par     - real vector                      : parameters, par=[] if no parameters
%        derivativeChoice   - string                : 'eqvar' or 'ind'
%
%    Outputs
%
%        tout    - real row vector                  : time at each integration step
%        z       - real matrix                      :  z(:,i) : flow at tout(i)
%        dz      - real matrix                      : dz(:,i) : linearized flow at tout(i)
%
%-------------------------------------------------------------------------------------------

X0  = [z0; dz0(:)];
n   = length(z0)/2;
k   = length(dz0)/(2*n);

odefun          = @(t,X) expdhvrhs(t, X, par, derivativeChoice, n, k);
[tout, Xout]    = ode45(odefun, tspan, X0, options);
tout            = tout';
z               = Xout(:,1:2*n)';
dz              = Xout(:,2*n+1:end)';

return

% On integre z et delta z en meme temps !
function rhs = expdhvrhs(t, X, par, derivativeChoice, n, k)

    rhs = zeros(length(X),1);
    z   = X(1:2*n);
    rhs(1:2*n) = hvfun(t,z,par);
    dhv = dhvfun(t, z, reshape(X(2*n+1:end),[2*n,k]), par, derivativeChoice);
    rhs(2*n+1:end) = dhv(:);

return
