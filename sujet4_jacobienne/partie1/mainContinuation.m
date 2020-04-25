% Auteur : Olivier Cots, INP-ENSEEIHT & IRIT
% Date   : Novembre 2016
%
% Etude du problème de contrôle optimal :
% min J(u) = int_{t0}^(tf} |u(t)| - epsilon * ( ln(|u(t)|) + ln(1-|u(t)|) ) dt
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
set(0,  'defaultaxesfontsize'   ,  14     , ...
'DefaultTextVerticalAlignment'  , 'bottom', ...
'DefaultTextHorizontalAlignment', 'left'  , ...
'DefaultTextFontSize'           ,  12     , ...
'DefaultFigureWindowStyle','docked');

figRes = figure('units','normalized'); %Visualisation
nbFigs = {2,2};

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
epsilon = 0; %Pas utile ici puisque mis à jour ensuite
iepsi = 5;


par = [t0; tf; x0; xf; epsilon];

n   = length(x0); % Dimension de l'état

tspan = [t0 tf];

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
optionsNLE              = optimoptions('fsolve','display','iter');%,'PlotFcn',@plotOptim);
options.odeset          = optionsODE;
options.optimoptions    = optionsNLE;

%-------------------------------------------------------------------------------------------------------------
% Exercice 2.1
%-------------------------------------------------------------------------------------------------------------
% Visualisation de la fonction de tir et de la solution, pour diff\'erentes valeurs de epsilon
%
disp(tirets);
disp('Exercice 2.1 : Visualisation de la fonction de tir et de la solution');

epsilons = [1 0.1 0.01 0.001 1e-6]; %On affiche pour plusieurs valeur de epsilon
colors   = {'b','y','k','g','r'};
Neps     = length(epsilons);

Npoints = 200;
yspan = linspace(0.0,1.5,Npoints); % p0 in [-2, 2]

leg   = {};
ysols = [];
ssols = [];
for i=1:Neps
    epsilon     = epsilons(i);
    par(iepsi)  = epsilon;

    % Calcul et affichage de la fonction de tir
    s     = [];
    for j=1:Npoints
        s  = [s sfun(yspan(j),optionsODE,par)];
    end;
    subplot(nbFigs{:},4); hold on; grid on; plot(yspan,s,colors{i},'LineWidth',LW);

    % Calcul de la solution
    %
    % Essayer y0 = 0.37, 0.38 et 0.39. Valeur par defaut : y0 = 0.2715
    %
    y0 = 0.2715;

    [ ysol, ssol, nfev, niter, flag] = ssolve(y0,options,par); ysols(i) = ysol; ssols(i) = ssol;

    fprintf(['Epsilon courant               : epsi  = %g\n'], epsilon);
    fprintf(['Solution trouvée              : p0    = %g\n'], ysol);
    fprintf(['Valeur de la fonction de tir  : S(p0) = %g\n'], ssol);
    fprintf(['Nombre évaluations de S       : nfev  = %i\n'], nfev);
    fprintf(['Nombre itérations             : niter = %i\n'], niter);

    % On affiche la solution
    z0          = [x0; ysol];
    [ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

    u = [];
    for j=1:length(tout)
        u = [u ; control(tout(j), z(:,j), par)];
    end;

    subplot(nbFigs{:},1); hold on; plot(tout, z(1,:), colors{i}, 'LineWidth', LW);
    subplot(nbFigs{:},2); hold on; plot(tout, z(2,:), colors{i}, 'LineWidth', LW);
    subplot(nbFigs{:},3); hold on; plot(tout,      u, colors{i}, 'LineWidth', LW);

    % Pour la legende
    leg{i} = sprintf(['epsilon = %g'], epsilon);

    % On affiche
    drawnow;
end;

% Les légendes
subplot(nbFigs{:},1); legend(leg{:},'Location','NorthWest');
subplot(nbFigs{:},2); legend(leg{:},'Location','NorthWest');
subplot(nbFigs{:},3); legend(leg{:},'Location','NorthWest');
subplot(nbFigs{:},4); legend(leg{:},'Location','SouthEast');
subplot(nbFigs{:},1); daxes(tf,xf,'k--');

% Les labels
subplot(nbFigs{:},1); xlabel('t'); ylabel('$x_\varepsilon$','Interpreter','LaTex');
subplot(nbFigs{:},2); xlabel('t'); ylabel('$p_\varepsilon$','Interpreter','LaTex');
subplot(nbFigs{:},3); xlabel('t'); ylabel('$u_\varepsilon$','Interpreter','LaTex');
subplot(nbFigs{:},4); xlabel('$p_0$','Interpreter','LaTex'); ylabel('$S_\varepsilon$','Interpreter','LaTex');

subplot(nbFigs{:},4); hold on; ax = axis; plot([ax(1) ax(2)], [0.0 0.0], 'k--', 'LineWidth',1.0);
for i=1:Neps
    ysol = ysols(i);
    ssol = ssols(i);
    subplot(nbFigs{:},4); hold on; plot(ysol, ssol, [colors{i} 'd'], 'MarkerFaceColor', colors{i}); % S(ysol)
    ax = axis; plot([ysol ysol], [ax(3) ax(4)], 'k--', 'LineWidth',1.0);
end;

return