function u = control(t, z, par, iarc)
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
%        u = control(t, z, par, iarc)
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
%        u   -  real vector : control
%
%-------------------------------------------------------------------------------------------
n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);

%% A REMPLACER
u = 1.0
if iarc == 1
    u = 0.0
end

