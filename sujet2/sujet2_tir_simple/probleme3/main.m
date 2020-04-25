
% Auteur : Olivier Cots, INP-ENSEEIHT & IRIT
% Date   : Novembre 2016
%
% Etude du problème de contrôle optimal :
% min J(u) = 1/2 * int_{t0}^(tf} u(t)^2 dt
% dot{x1}(t) = x2()
% dot{x2}(t) = u(t), u(t) in [-umax, umax]
% x(0) = x0
% x(tf) = xf
%
% t0 = 0, tf = 1, x0 = (-1,0), xf = (0,0).
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
nbFigs = {2,3};

% On choisit un format long pour plus de détails sur les résultats numériques
%
format long;

% On ajoute dans le path, le répertoire lib contenant les routines à implémenter par l'étudiant
%
addpath(['lib/']);
addpath(['ressources/']);

% Définition des paramètres : par = [t0, tf, x0, xf, umax]
% Et du tspan
%
t0  =  0.0;
tf = 1.0;
x0 = [-1.0; 0.0];
xf = [ 0.0; 0.0];
umax= 7.0;
par = [t0; tf; x0(:); xf(:); umax];

n = length(x0); % Dimension de l'état

tspan = [t0 tf];

% Définitions des options
%
optionsODE = odeset;
optionsNLE = optimoptions('fsolve','display','iter');%,'PlotFcn',@plotOptim);
options.odeset = optionsODE;
options.optimoptions = optionsNLE;

%-------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------
% Affichage d'intégrations numériques
%
% p0_1
p0_1 = [1.0; 1.0];
z0 = [x0; p0_1];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z( 1,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('x_1');
subplot(nbFigs{:},2); hold on; plot(tout, z( 2,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('x_2');
subplot(nbFigs{:},3); hold on; plot(tout, u, 'b', 'LineWidth', LW); xlabel('t'); ylabel('u');
subplot(nbFigs{:},4); hold on; plot(tout, z(n+1,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('p_1');
subplot(nbFigs{:},5); hold on; plot(tout, z(n+2,:), 'b', 'LineWidth', LW); xlabel('t'); ylabel('p_2');

% p0_2
p0_2 = [10.0; 10.0];
z0 = [x0; p0_2];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z( 1,:), 'r', 'LineWidth', LW); xlabel('t'); ylabel('x_1');
subplot(nbFigs{:},2); hold on; plot(tout, z( 2,:), 'r', 'LineWidth', LW); xlabel('t'); ylabel('x_2');
subplot(nbFigs{:},3); hold on; plot(tout, u, 'r', 'LineWidth', LW); xlabel('t'); ylabel('u');
subplot(nbFigs{:},4); hold on; plot(tout, z(n+1,:), 'r', 'LineWidth', LW); xlabel('t'); ylabel('p_1');
subplot(nbFigs{:},5); hold on; plot(tout, z(n+2,:), 'r', 'LineWidth', LW); xlabel('t'); ylabel('p_2');

% Les légendes et autres
subplot(nbFigs{:},1); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2))},'Location','NorthWest');
subplot(nbFigs{:},2); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2))},'Location','NorthWest');
subplot(nbFigs{:},3); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2))},'Location','NorthEast');
subplot(nbFigs{:},4); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2))},'Location','NorthWest');
subplot(nbFigs{:},5); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2))},'Location','NorthEast');

%-------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------
% Affichage de la norme de la fonction de tir
%
NbPoints = 30;
X1 = linspace(5,15,NbPoints);
X2 = linspace(0,10,NbPoints);
[P01, P02] = meshgrid(X1,X2);
s = [];
fprintf(['Calcul de la fonction de tir sur la grille %g x %g\n'], NbPoints, NbPoints);
for i=1:NbPoints
for j=1:NbPoints
s(i,j) = norm(sfun([P01(i,j); P02(i,j)], optionsODE, par));
fprintf(1, '\b\b\b\b\b\b%3.1f%% ', 100*((i-1)*NbPoints+j)/(NbPoints^2));
end;
end;
fprintf(1,'\n');

subplot(nbFigs{:},6); hold on; grid on; surf(P01, P02, s,'LineWidth',LW); xlabel('p_0(1)'); ylabel('p_0(2)'); zlabel('||S||');
az = -45; el = 45; view(az, el);

% Pause
disp('Appuyer sur une touche pour continuer...');
pause

%-------------------------------------------------------------------------------------------------------------
% Exercice 1.6
%-------------------------------------------------------------------------------------------------------------
% Recherche d'une extrémale vérifiant les conditions nécessaires d'optimalité,
% ie d'une BC-extrémale. Résolution de S(y) = 0 par la méthode de tir (Newton).
%
disp(tirets);
disp('Exercice 1.6 : problème en dimension 2 vec contraintes sur le contrôle');

%% Valeur initiale pour la méthode de tir
y0 = [5; 5];

[ ysol, ssol, nfev, niter, flag] = ssolve(y0,options,par);

fprintf(['Solution trouvée : p0 = (%g, %g)\n'], ysol(1), ysol(2));
fprintf(['Valeur de la fonction de tir : S(p0) = (%g, %g)\n'], ssol(1), ssol(2));
fprintf(['Nombre évaluations de S : nfev = %i\n'], nfev);
fprintf(['Nombre itérations : niter = %i\n\n'], niter);

% On affiche la solution
figure(figRes);
subplot(nbFigs{:},6); hold on; daxes3(ysol(1),ysol(2),norm(ssol),'k--');

z0 = [x0; ysol];
[ tout, z ] = exphvfun(tspan, z0, optionsODE, par);

% calcul du contrôle
u = [];
for i=1:length(tout)
u = [u ; control(tout(i), z(:,i), par)];
end;

subplot(nbFigs{:},1); hold on; plot(tout, z( 1,:), 'm', 'LineWidth', LW); xlabel('t'); ylabel('x_1');
subplot(nbFigs{:},2); hold on; plot(tout, z( 2,:), 'm', 'LineWidth', LW); xlabel('t'); ylabel('x_2');
subplot(nbFigs{:},3); hold on; plot(tout, u, 'm', 'LineWidth', LW); xlabel('t'); ylabel('u');
subplot(nbFigs{:},4); hold on; plot(tout, z(n+1,:), 'm', 'LineWidth', LW); xlabel('t'); ylabel('p_1');
subplot(nbFigs{:},5); hold on; plot(tout, z(n+2,:), 'm', 'LineWidth', LW); xlabel('t'); ylabel('p_2');

% Les légendes et autres
subplot(nbFigs{:},1); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2)),sprintf(['p_0 = solution'])},'Location','NorthWest');
subplot(nbFigs{:},2); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2)),sprintf(['p_0 = solution'])},'Location','NorthWest');
subplot(nbFigs{:},3); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2)),sprintf(['p_0 = solution'])},'Location','NorthEast');
subplot(nbFigs{:},4); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2)),sprintf(['p_0 = solution'])},'Location','NorthWest');
subplot(nbFigs{:},5); legend({sprintf(['p_0 = (%g,%g)'],p0_1(1),p0_1(2)),sprintf(['p_0 = (%g,%g)'],p0_2(1),p0_2(2)),sprintf(['p_0 = solution'])},'Location','NorthEast');

% Cible finale
subplot(nbFigs{:},1); daxes(1,0,'k--');
subplot(nbFigs{:},2); daxes(1,0,'k--');
