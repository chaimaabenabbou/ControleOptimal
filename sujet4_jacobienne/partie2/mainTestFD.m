% Auteur : Olivier Cots, INP-ENSEEIHT & IRIT
% Date   : Novembre 2016
%
% Test de l'algorithme de différences finies avant
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
nbFigs = {1,1};

% On choisit un format long pour plus de détails sur les résultats numériques
%
format long;

% On ajoute dans le path, le répertoire lib contenant les routines à implémenter par l'étudiant
%
addpath(['lib/']);
addpath(['ressources/']);

%-------------------------------------------------------------------------------------------------------------
% Exercice 2.2
%-------------------------------------------------------------------------------------------------------------
% Test de l'algorithme de différences finies avant sur la fonction cosinus
%
disp(tirets);
disp('Exercice 2.2 : Test de l''algorithme de différences finies avant sur la fonction cosinus');

fun = @(x) -cos(x);
jac = @(x)  sin(x);

M1  = 1.0;
Lf  = 1.0;
u   = eps;
phi = @(t)  M1*t/2+2*u*Lf./t;

NbPoints = 100;
xspan    = linspace(0,2*pi,NbPoints);

I        = eye(1);

erreurs  = [];
ordres   = [1:16];
for i = 1:length(ordres)

    t = 10^(-ordres(i));

    jacFinite = [];
    for j=1:NbPoints
        jacFinite(j) = finiteDiff(fun,xspan(j),I,t);
    end;

    erreurs(i) = mean(abs(jacFinite - jac(xspan)));

end;

subplot(nbFigs{:},1);
loglog(10.^(-ordres), erreurs, 'b', 'LineWidth', LW); hold on; grid on; xlabel('t'); ylabel('erreurs en moyenne');
loglog(10.^(-ordres), phi(10.^(-ordres)), 'r', 'LineWidth', LW);
title('Erreur en moyenne des differences finies avants sur la fonction cosinus');
legend({'Numerique','Theorique'},'Location','NorthEast');
