function [ J ] = expcost(tspan, z0, options, par, iarc)
%-------------------------------------------------------------------------------------------
%
%    expcost
%
%    Description
%
%        Computes the cost.
%
%-------------------------------------------------------------------------------------------
%
%    Matlab Usage
%
%        [tout, z, flag] = expcost(tspan, z0, options, iepsi)
%
%    Inputs
%
%        tspan   - real row vector of dimension 2   : tspan = [t0 tf]
%        z0      - real vector                      : initial flow
%        options - struct vector                    : odeset options
%        par     - real vector                      : parameters, par=[] if no parameters
%        iarc    -  integer                         : index of current arc
%
%    Outputs
%
%        J       - real                             : cost
%
%-------------------------------------------------------------------------------------------

odefun      = @(t,zc) rhsfun(t,zc,par,iarc);
[tout, zc]  = ode45(odefun, tspan, [z0; 0.0], options);
J           = zc(end);

return

function rhs = rhsfun(t, zc, par, iarc)

n   = (length(zc)-1)/2;
z   = zc(1:2*n);
u   = control(t,z,par,iarc);
rhs = zeros(length(zc),1);
rhs(1:2*n) = hvfun(t,z,par,iarc);
au         = abs(u);
rhs(2*n+1) = au; % |u(t)|

return
