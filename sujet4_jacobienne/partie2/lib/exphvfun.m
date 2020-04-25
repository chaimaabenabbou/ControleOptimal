function [ tout, z ]  = exphvfun(tspan, z0, options, par)
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
%        [tout, z, flag] = exphvfun(tspan, z0, options, par)
%
%    Inputs
%
%        tspan   - real row vector of dimension 2   : tspan = [t0 tf]
%        z0      - real vector                      : initial flow
%        options - struct vector                    : odeset options
%        par     - real vector                      : parameters, par=[] if no parameters
%
%    Outputs
%
%        tout    - real row vector                  : time at each integration step
%        z       - real matrix                      : z(:,i) : flow at tout(i)
%
%-------------------------------------------------------------------------------------------

%% A REMPLACER
tout = linspace(tspan(1), tspan(2), 100);
z    = z0 * ones(1,length(tout));
