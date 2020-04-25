function dhv = dhvfun(t, z, dz, par, derivativeChoice)
%-------------------------------------------------------------------------------------------
%
%    dhvfun
%
%    Description
%
%        Computes the linearized of the Hamiltonian vector field associated to H.
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        dhv = dhvfun(t, z, dz, par, derivativeChoice)
%
%    Inputs
%
%        t                  -  real        : time
%        z                  -  real vector : state and costate
%        dz                 -  real matrix : variation of z
%        par                -  real vector : parameters, par=[] if no parameters
%        derivativeChoice   -  string      : 'eqvar' or 'ind'
%
%    Outputs
%
%        dhv  -  real vector : dhv = dhv(z) . dz
%
%-------------------------------------------------------------------------------------------

switch derivativeChoice

    case 'ind' % Diff\'erence finie sur le syst\`eme hamiltonien

        dhv  = zeros(size(dz));

    case 'eqvar'

        dhv  = zeros(size(dz));

    otherwise
        error('derivativeChoice ne peut prendre les valeurs : ''eqvar'' ou ''ind''');
end;
