function [ tout, z ]  = exphvfun(tspan, z0, options, par, iarc)
%-------------------------------------------------------------------------------------------
%
%    exphvfun
%
%    Description
%
%        Computes the chronological exponential of the Hamiltonian vector field hv
%        defined by h.
%
%-------------------------------------------------------------------------------------------
%
%    Matlab Usage
%
%        [tout, z, flag] = exphvfun(tspan, z0, options, par, iarc)
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
%        tout    - real row vector                  : time at each integration step
%        z       - real matrix                      : z(:,i) : flow at tout(i)
%
%-------------------------------------------------------------------------------------------

%% A REMPLACER
[tout, z] = ode45(@(t,z) hvfun(t, z, par, iarc), tspan, z0,options);
tout = tout';
z=z';
