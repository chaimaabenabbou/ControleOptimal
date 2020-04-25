function s  = sfun(y,options,par)
%-------------------------------------------------------------------------------------------
%
%    sfun
%
%    Description
%
%        Computes the shooting function
%
%-------------------------------------------------------------------------------------------
%
%    Matlab
%
%        s = sfun(y, options, par)
%
%    Inputs
%
%        y       - real vector  : shooting variable
%        options - struct       : odeset options
%        par     - real vector  : parameters, par=[] if no parameters
%
%    Outputs
%
%        s       - real vector, shooting value
%
%-------------------------------------------------------------------------------------------

% par = [t0; tf; x0; xf, epsilon]
%
t0 = par(1);
tf = par(2);
x0 = par(3);
xf = par(4);

p0 = y;
z0 = [x0;p0];

%% A REMPLACER
s = zeros(length(y),1);

