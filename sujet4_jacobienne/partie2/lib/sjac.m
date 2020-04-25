function jac  = sjac(y,options,par,derivativeChoice)
%-------------------------------------------------------------------------------------------
%
%    sjac
%
%    Description
%
%        Computes the Jacobian of the shooting function
%
%-------------------------------------------------------------------------------------------
%
%    Matlab
%
%        s = sjac(y, options, par, derivativeChoice)
%
%    Inputs
%
%        y                  - real vector  : shooting variable
%        options            - struct       : odeset options
%        par                - real vector  : parameters, par=[] if no parameters
%        derivativeChoice   - String       : must be 'finite', 'ind' or 'eqvar'
%           if 'finite' then sjac      is computed by finite differences
%           if 'ind'    then expdhvfun is computed by intern finite differences (BOCK)
%           if 'eqvar'  then expdhvfun is computed by variational equations
%
%    Outputs
%
%        jac     - real matrix, Jacobian of the shooting value
%
%-------------------------------------------------------------------------------------------

if(nargin==3)
    derivativeChoice = 'eqvar';
end;

switch derivativeChoice

    case 'finite' %Calcul de la jacobienne de la fonction de tir par diff\'erences finies avants

        abstol = options.AbsTol;
        reltol = options.RelTol;
        if(isempty(abstol) & isempty(reltol))
            t = sqrt(eps);          % On copie le comportement par d\'efaut de la fonction fsolve
                                    % de matlab qui utilise un pas de l'ordre de la racine carr\'ee
                                    % de l'epsilon machine pour l'approximation par diff\'erences
                                    % finies.
            %t = sqrt(1e-3);
        else
            t = sqrt(max(abstol,reltol)); % on choisit un pas de l'ordre de la racine carr\'ee de
                                    % l'erreur num\'erique du calcul de la jacobienne, qui ici
                                    % vient de l'int\'egration num\'erique du syst\`eme hamiltonien
        end;

        %% A REMPLACER
        jac = zeros(length(y), length(y));

    case {'ind', 'eqvar'}

        % par = [t0; tf; x0; xf, epsilon]
        %

        %% A REMPLACER
        jac = zeros(length(y), length(y));

    otherwise
        error('Choix ne peut prendre les valeurs : ''finite'', ''eqvar'' ou ''ind''');
end
