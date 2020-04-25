% Auteur : Olivier Cots, INP-ENSEEIHT & IRIT
% Date   : Novembre 2016
%
% Etude du problème de contrôle optimal :
% min J(u) = 1/2 * int_{t0}^(tf} u(t)^2 dt
% dot{x}(t) = - x(t) + u(t), u(t) in R
% x(0) = x0
% x(tf) = xf
%
% t0 = 0, tf = 1, x0 = -1, xf = 0.
%

% On réinitialise l'environnement
%
clear all;
close all;
path(pathdef);

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
nbFigs = {2,2};

% On choisit un format long pour plus de détails sur les résultats numériques
%
format long;

% On ajoute dans le path, le répertoire lib contenant les routines à implémenter par l'étudiant
%
addpath(['lib/']);
addpath(['ressources/']);

% Définition des paramètres : par = [t0, tf, x0, xf]
% Et du tspan
%
t0  =  0.0;
tf  =  1.0;
x0  = -1.0;
xf  =  0.0;
par = [t0; tf; x0; xf];

tspan = [t0 tf];

% Définitions des options
%
optionsODE              = odeset;
optionsNLE              = optimoptions('fsolve','display','iter');%,'PlotFcn',@plotOptim);
options.odeset          = optionsODE;
options.optimoptions    = optionsNLE;

% Solution
%
p0sol = 2.0*(xf*exp(tf)-x0)/(exp(2.0*tf)-1.0);

%-------------------------------------------------------------------------------------------------------------
% Exercice 1.1
%-------------------------------------------------------------------------------------------------------------
% Vérification de l'implémentation de hvfun
%
disp(tirets);
disp('Exercice 1.1 : hvfun');

p0  = 0.0;
z0  = [x0; p0];
hv  = hvfun(t0, z0, par);
fprintf(['hvfun(0.0, [-1.0; 0.0], par) = (%g,%g)\n'], hv(1),hv(2));

%-------------------------------------------------------------------------------------------------------------
% Exercice 1.2
%-------------------------------------------------------------------------------------------------------------
% Vérification de l'intégration numérique
%
disp(tirets);
disp('Exercice 1.2 : exphvfun');

% p0 = 0.0
p0_1        = 0.0;
z0          = [x0; p0_1];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
    u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z(1,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('x');
subplot(nbFigs{:},2); hold on; plot(tout, z(2,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('p');
subplot(nbFigs{:},3); hold on; plot(tout,      u, 'b', 'LineWidth', LW); xlabel('t'); ylabel('u');

% p0 = 0.4
p0_2        = 0.4;
z0          = [x0; p0_2];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
    u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z(1,:), 'r', 'LineWidth', LW);
subplot(nbFigs{:},2); hold on; plot(tout, z(2,:), 'r', 'LineWidth', LW);
subplot(nbFigs{:},3); hold on; plot(tout,      u, 'r', 'LineWidth', LW);

% Les légendes et autres
subplot(nbFigs{:},1); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2)},'Location','NorthWest');
subplot(nbFigs{:},2); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2)},'Location','NorthWest');
subplot(nbFigs{:},3); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2)},'Location','NorthWest');

%-------------------------------------------------------------------------------------------------------------
% Exercice 1.3
%-------------------------------------------------------------------------------------------------------------
% Vérification de la fonction de tir
%
disp(tirets);
disp('Exercice 1.3 : sfun');

yspan = linspace(-3,3,100); % p0 in [-3, 3]
s     = [];
for i=1:length(yspan)
    s  = [s sfun(yspan(i),optionsODE,par)];
end;

subplot(nbFigs{:},4); hold on; plot(yspan,s,'b','LineWidth',LW); xlabel('p_0'); ylabel('S'); daxes(p0sol,0,'k--');

% Pause
disp('Appuyer sur une touche pour continuer...');
pause

%-------------------------------------------------------------------------------------------------------------
% Exercice 1.4
%-------------------------------------------------------------------------------------------------------------
% Recherche d'une extrémale vérifiant les conditions nécessaires d'optimalité,
% ie d'une BC-extrémale. Résolution de S(y) = 0 par la méthode de tir (Newton).
%
disp(tirets);
disp('Exercice 1.4 : ssolve ou la méthode de tir');

y0 = 0.5;
[ ysol, ssol, nfev, niter, flag] = ssolve(y0,options,par);

fprintf(['Solution trouvée              : p0    = %g\n'], ysol);
fprintf(['Ecart à la solution           : e     = %g\n'], ysol-p0sol);
fprintf(['Valeur de la fonction de tir  : S(p0) = %g\n'], ssol);
fprintf(['Nombre évaluations de S       : nfev  = %i\n'], nfev);
fprintf(['Nombre itérations             : niter = %i\n\n'], niter);

% On affiche la solution
figure(figRes);
subplot(nbFigs{:},4); hold on; plot(ysol, ssol,'kd', 'MarkerFaceColor',[1 0 1]); % S(ysol)
axis([min(yspan) max(yspan) min(s)-0.01 max(s)+0.01]) % S(ysol)

z0          = [x0; ysol];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
    u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z(1,:), 'm', 'LineWidth', LW);
subplot(nbFigs{:},2); hold on; plot(tout, z(2,:), 'm', 'LineWidth', LW);
subplot(nbFigs{:},3); hold on; plot(tout,      u, 'm', 'LineWidth', LW);

% Les légendes et autres
subplot(nbFigs{:},1); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2),sprintf(['p_0 = %g'],ysol)},'Location','NorthWest');
subplot(nbFigs{:},2); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2),sprintf(['p_0 = %g'],ysol)},'Location','NorthWest');
subplot(nbFigs{:},3); legend({sprintf(['p_0 = %g'],p0_1),sprintf(['p_0 = %g'],p0_2),sprintf(['p_0 = %g'],ysol)},'Location','NorthWest');
subplot(nbFigs{:},1); daxes(1,0,'k--');
