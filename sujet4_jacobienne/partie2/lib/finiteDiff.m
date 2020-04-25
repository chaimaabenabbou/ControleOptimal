function dfun = finiteDiff(fun,x,h,t)
%-------------------------------------------------------------------------------------------
%
%    finiteDiff
%
%    Description
%
%        Computes an approximation of the directional derivative of fun at x along h,
%        by forward finite differences.
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        dfun = finiteDiff(fun,x,h,t)
%
%    Inputs
%
%        fun  -  the function
%        x    -  real vector : the current point
%        h    -  real matrix : the directions along which the differential is computed,
%                               h = [h_1 h_2 ... h_m]
%        t    -  real scalar : the absolute step for approximation
%
%    Outputs
%
%        dfun  -  real matrix : approximation of dfun(x) . h
%
%                               dfun(:,i) = (fun(x+t*h_i)-fun(x))/t, h = [h_1 h_2 ... h_m]
%
%-------------------------------------------------------------------------------------------
m = length(h);

for i = 1:m
   h_i = h(i);
   dfun(:,i)=(fun(x+t*h_i) - fun(x))/t;    
end

