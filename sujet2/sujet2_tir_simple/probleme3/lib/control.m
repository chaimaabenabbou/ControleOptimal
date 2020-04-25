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
% par = [t0; tf; x0; xf; umax]
%
n  = length(z)/2;
x  = z(1:n);
p  = z(n+1:2*n);

p2 = p((length(p)/2)+1:end);
umax = par(end);
%% A REMPLACER
for k=1:length(p2)
    if abs(p2(k))<=umax
        u(k)=p2(k);
    elseif p2(k)>umax
        u(k)=umax;
    else
        u(k)=-umax;
    end
end






