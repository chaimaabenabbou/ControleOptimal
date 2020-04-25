%
% ~gergaud/ENS/ontrole/TP-19-20/Codes_Etudiant/sujet1_RKE/IVP_pb1_conso.m
%
% Auteurs:  Joseph GERGAUD
% Date:     février 2020
% Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
%           2, rue Camichel 31071 Toulouse FRANCE
% Email:    gergaud@enseeiht.fr
%***************************************************************************
%
% Intégration de l'equation differentiel du problème 1 avec minimisation de
% la consommation.
%
%
% Initialisations
% ---------------
clear all; 
close all;
restoredefaultpath;

addpath('./lib');

format long

hvfun = @hvfun_pb1_conso;
exphvfun = @exphvfun_pb1_conso;
x0 = 0; p0 = exp(-1);    % ou p0 = exp(-1.1)

z0=[x0 ; p0];
t0=0;
tf=2.1;
N0=[120:60:1080 1200:600:10800 ]; % Valeurs de N pour les courbes d'ordre
%
% Solutions z_1 et z_2 et plan de phase pour RKE
% ----------------------------------------------
N = 10;
figure(1)
disp('Euler')
disp('-----')
[T,Z]=ode_euler(hvfun,[t0 tf],z0,N);
disp('[T Z]=');
disp([T Z])
plot_sol(T,Z)

%%%%%%%%%%%%%%%%%
disp('Runge')
disp('-----')
[T,X] = ode_runge(hvfun,[t0 tf],z0,N);
disp('[T X]=');
disp([T X])
plot_sol(T,X)

%%%%%%%%%%%%%%%%
disp('Heun')
disp('-----')
[T,Xheun] = ode_heun(hvfun,[t0 tf],z0,N);
disp('[T Xheun]=');
disp([T Xheun])
plot_sol(T,Xheun)


%%%%%%%%%%%%%%%%%%
disp('Rk4')
disp('-----')
[T,Xrk4] = ode_rk4(hvfun,[t0 tf],z0,N);
disp('[T Xrk4]=');
disp([T Xrk4])
plot_sol(T,Xrk4)

legend('Euler','Runge', 'Heun', 'RK4','ode45');

%% ORDRE
%print('fig_pb1_cons','-depsc')

%Courbes d'ordre
%---------------
N=N0;
%erreur (ordre)
zf = exphvfun_pb1_conso(tf,z0);
C1 = zeros(length(N),2);
C2 = C1;
C3 = C1;
C4 = C1;
nfe= N;

%Euler
for i = 1:length(N)
    [T,Z]=ode_euler(hvfun,[t0 tf],z0,N(i)/1);
    C1(i,1) = log10(abs( Z(end,1)' - zf(1))) ;
    C1(i,2) = log10(abs( Z(end,2)' - zf(2))) ;
end

%Runge
for i = 1:length(N)
    [T,Z]=ode_runge(hvfun,[t0 tf],z0,N(i)/2);
    C2(i,1) = log10(abs( Z(end,1)' - zf(1))) ;
    C2(i,2) = log10(abs( Z(end,2)' - zf(2))) ;
end

%Heun
for i = 1:length(N)
    [T,Z]=ode_heun(hvfun,[t0 tf],z0,N(i)/3);
    C3(i,1) = log10(abs( Z(end,1)' - zf(1))) ;
    C3(i,2) = log10(abs( Z(end,2)' - zf(2))) ;
end

%RK4
for i = 1:length(N)
    [T,Z]=ode_rk4(hvfun,[t0 tf],z0,N(i)/4);
    C4(i,1) = log10(abs( Z(end,1)' - zf(1))) ;
    C4(i,2) = log10(abs( Z(end,2)' - zf(2))) ;
end

figure(2)
subplot(1,2,1)
plot(log10(nfe),C1(:,1))
hold on 
plot(log10(nfe),C2(:,1))
hold on
plot(log10(nfe),C3(:,1))
hold on
plot(log10(nfe),C4(:,1))

xlabel('log10(fe)') 
ylabel('log10(erreur pour z1)') 

subplot(1,2,2)
plot(log10(nfe),C1(:,2))
hold on 
plot(log10(nfe),C2(:,2))
hold on
plot(log10(nfe),C3(:,2))
hold on
plot(log10(nfe),C4(:,2))

xlabel('log10(fe)') 
ylabel('log10(erreur pour z2)') 

legend('Euler','Runge', 'Heun', 'RK4');


% -----------------------------------------------------------------------
% deuxieme membre de l'edo
% ------------------------
%
function zpoint = hvfun_pb1_conso(t,z)
%
% A remplacer
 
    if abs(z(2))<1 
        u=0;
    else
        u = sign(z(2));
    end
    zpoint = [-z(1) + u ; z(2)];
end

%
% Solution exacte de l'edo
% ------------------------
function zf = exphvfun_pb1_conso(t,z0)
%
% A remplacer
    zf2=z0(2)*exp(t);
    zf1=((z0(1)-exp(1))*exp(-t(length(t))))+1;
    zf=[zf1,zf2(end)];
end
%
% Fonction auxiliaires
% --------------------
function plot_sol(T,Z);
% plot les soltions
% T = temps;
% Y = solutions
% plot les soltions
% T = temps;
% Y = solutions
    dim = size(T);
    if dim(1) == 1
      T = T(:);
      Z = Z';
    end  

    subplot(2,2,1)
    hold on
    plot(T,Z(:,1))
    xlabel('t')
    ylabel('z_1(t)')
    subplot(2,2,2)
    hold on
    plot(T,Z(:,2))
    xlabel('t')
    ylabel('z_2(t)')
    subplot(2,2,3:4)
    hold on
    plot(Z(:,1),Z(:,2))
    xlabel('z_1(t)')
    ylabel('z_2(t)')
end