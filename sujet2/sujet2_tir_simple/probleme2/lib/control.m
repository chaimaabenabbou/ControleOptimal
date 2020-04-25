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
n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);

%% A REMPLACER
%u = 0.0;
if abs(p) <= 1 
    u = p;
else
    u =  sign(p);
    
end