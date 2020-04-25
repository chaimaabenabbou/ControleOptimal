% Auteur : Olivier Cots, INP-ENSEEIHT & IRIT
% Date   : mai 2017
%
% Résolution d'un problème bang-bang:
% min J(u) = int_{t0}^(tf} |u(t)| dt
% dot{x}(t) = -x(t) + u(t), u(t) in [-1, 1]
% x(0) = x0
% x(tf) = xf
%
% t0 = 0, tf = 2, x0 = 0, xf = 0.5.
%

% On réinitialise l'environnement
%
clear all;
close all;
path(pathdef);
clc;

% Des paramètres pour l'affichage
%
tirets  = ['------------------------------------------'];
LW      = 1.5;
set(0,  'defaultaxesfontsize'   ,  12     , ...
'DefaultTextVerticalAlignment'  , 'bottom', ...
'DefaultTextHorizontalAlignment', 'left'  , ...
'DefaultTextFontSize'           ,  12     , ...
'DefaultFigureWindowStyle','docked');

figRes = figure('units','normalized'); %Visualisation
nbFigs = {1,3};

% On choisit un format long pour plus de détails sur les résultats numériques
%
format long;

% On ajoute dans le path, le répertoire lib contenant les routines à implémenter par l'étudiant
%
addpath(['lib/']);
addpath(['ressources/']);

% Définition des paramètres : par = [t0, tf, x0, xf, epsilon]
% Et du tspan
%
t0  = 0.0;
tf  = 2.0;
x0  = 0.0;
xf  = 0.5;

n   = length(x0); % Dimension de l'état

% Définitions des options
%
optionParDefaut = 1;
if(optionParDefaut)
    Tol                 = 1e-3;
    optionsODE          = odeset; %('AbsTol',1e-6,'RelTol',1e-3);
else
    Tol                 = 1e-6;
    optionsODE          = odeset('AbsTol',Tol,'RelTol',Tol);
end;
optionsNLE              = optimoptions('fsolve','display','iter','FiniteDifferenceStepSize',sqrt(Tol));
options.odeset          = optionsODE;
options.optimoptions    = optionsNLE;
%-------------------------------------------------------------------------------------------------------------
% Exercice 4.1
%-------------------------------------------------------------------------------------------------------------
% Tir multiple et visualisation des résultats
%
disp(tirets);
disp('Exercice 4.1 : Tir multiple et visualisation des résultats.');


%--------------------------------
p0  = 0.270671;
t1  = 1.30685;
z1  = [0;1];
%--------------------------------

y0  = [p0; t1; z1];

par = [t0; tf; x0; xf];

[ ysol, ssol, nfev, niter, flag] = ssolve(y0,options,par);

p0 = ysol(1);
t1 = ysol(2);
z1 = ysol(3:4);

fprintf(['Solution trouvée              : p0      = %g\n'], p0);
fprintf(['                              : t1      = %g\n'], t1);
fprintf(['                              : z1      = (%g, %g)\n'], z1(1), z1(2));
fprintf(['Valeur de la fonction de tir  : S(ysol) = (%g, %g, %g, %g)\n'], ssol(1), ssol(2), ssol(3), ssol(4));
fprintf(['Nombre évaluations de S       : nfev    = %i\n'], nfev);
fprintf(['Nombre itérations             : niter   = %i\n'], niter);

% On affiche la solution
figure(figRes);

tout = [];
z    = [];
u    = [];
J    = 0;
for iarc = 1:2

    if(iarc == 1)
        tspan   = [t0 t1];
        z0      = [x0; p0];
    elseif(iarc==2)
        tspan   = [t1 tf];
        z0      = z1;
    end

    [ tout_aux, z_aux ] = exphvfun(tspan, z0, optionsODE, par, iarc);

    tout = [tout tout_aux];
    z    = [z    z_aux];

    % calcul du contrôle
    for i=1:length(tout_aux)
        u = [u ; control(tout_aux(i), z_aux(:,i), par, iarc)];
    end;

    % Calul du cout
    J = J + expcost(tspan, z0, optionsODE, par, iarc);

end;

fprintf(['Valeur du critère             : J       = %g\n\n'], J);

subplot(nbFigs{:},1); hold on; plot(tout, z(  1,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('x');
subplot(nbFigs{:},2); hold on; plot(tout, z(n+1,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('p');
subplot(nbFigs{:},3); hold on; plot(tout,        u, 'b', 'LineWidth', LW); xlabel('t'); ylabel('u');

% Cible finale
subplot(nbFigs{:},1); daxes(tf,xf,'k--');
subplot(nbFigs{:},3); daxes(0.0,1.0,'k--');

