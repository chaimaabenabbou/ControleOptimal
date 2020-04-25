function [ ysol, ssol, nfev, niter, flag] = ssolve(y0,options,par,derivativeChoice)
%-------------------------------------------------------------------------------------------
%
%    ssolve
%
%    Description
%
%        Interface of the Matlab non linear solver (fsolve) to solve the optimal
%        control problem described by the sfun subroutines.
%
%-------------------------------------------------------------------------------------------
%
%    Matlab Usage
%
%        [ysol,ssol,nfev,njev,flag] = ssolve(y0,options,par)
%
%    Inputs
%
%        y0      - real vector      : intial guess for shooting variable
%        options - struct vector    : options.odeset and options.optimoptions
%        par     - real vector      : parameters, par=[] if no parameters
%        derivativeChoice - String  : must be 'finite', 'ind' or 'eqvar'
%           if 'finite' then sjac      is computed by finite differences
%           if 'ind'    then expdhvfun is computed by intern finite differences (BOCK)
%           if 'eqvar'  then expdhvfun is computed by variational equations
%
%    Outputs
%
%        ysol    - real vector      : shooting variable solution
%        ssol    - real vector      : value of sfun at ysol
%        nfev    - integer          : number of evaluations of sfun
%        niter   - integer          : number of iterations
%        flag    - integer          : solver output (should be 1)
%
%-------------------------------------------------------------------------------------------

if(strcmp(options.optimoptions.Jacobian,'on')==1)
    if(nargin==3)
        derivativeChoice = 'eqvar';
    end;
    fun     = @(y) sfunjac(y,options.odeset,par,derivativeChoice);
else
    fun     = @(y) sfun(y,options.odeset,par);
end;

[ysol,ssol,flag,output,jacobian] = fsolve(fun,y0,options.optimoptions);

nfev    = output.funcCount;
niter   = output.iterations;

return

function [s,j] = sfunjac(y,options,par,derivativeChoice)

    s = sfun(y,options,par);
    if(nargout>1)
        j = sjac(y,options,par,derivativeChoice);
    end

return
