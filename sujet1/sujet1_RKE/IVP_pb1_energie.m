%
% ~gergaud/ENS/ontrole/TP-19-20/Codes_Etudiant/sujet1_RKE/IVP_pb1_energie.m
%
% Auteurs:  Joseph GERGAUD
% Date:     février 2020
% Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
%           2, rue Camichel 31071 Toulouse FRANCE
% Email:    gergaud@enseeiht.fr
%***************************************************************************
%
% Intégration de l'equation differentiel du problème 1 avec minimisation de
% l'énergie.
%
%
% Initialisations
% ---------------
clear all; close all;
restoredefaultpath;

addpath('./lib');

format long
%

hvfun = @hvfun_pb1_energie;                          % nom du deuxième membre de l'edo
exphvfun = @exphvfun_pb1_energie;                    % nom de la fonction solution de l'edo

x0 = 0; p0 = exp(-1);
z0=[x0 ; p0];
t0=0;
tf=2.1;
N0=[120:60:1080 1200:600:10800]; % Valeurs de N pour les courbes d'ordre
%
% Solutions x_1 et x_2 et plan de phase pour RKE
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
%print('fig_pb1_energie','-depsc')

%Courbes d'ordre
%---------------
N=N0;
%erreur (ordre)
zf = exphvfun_pb1_energie(tf,z0);
C1=zeros(length(N),2);
C2=C1;
C3=C1;
C4=C1;
C5=C1;

%Euler
nfe= N;
for i = 1:length(N)
    [T,Z]=ode_euler(hvfun,[t0 tf],z0,N(i)/1);
    C1(i,:) = log10(abs( Z(end,:)' - zf)) ;
end

%Runge
for i = 1:length(N)
    [T,Z]=ode_runge(hvfun,[t0 tf],z0,N(i)/2);
    C2(i,:) = log10(abs( Z(end,:)' - zf)) ;
end

%Heun
for i = 1:length(N)
    [T,Z]=ode_heun(hvfun,[t0 tf],z0,N(i)/3);
    C3(i,:) = log10(abs( Z(end,:)' - zf)) ;
end

% %RK4
for i = 1:length(N)
    [T,Z]=ode_rk4(hvfun,[t0 tf],z0,N(i)/4);
    C4(i,:) = log10(abs( Z(end,:)' - zf)) ;
end

% %%ode45
% RelTol = 1e-4; %[1e-4,1e-5,1e-6,1e-7,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-14,1e-15,1e-16];
% AbsTol = RelTol;   
% for i = 1:length(N)
%     options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
%     sol = ode45(hvfun,[t0,tf],z0,options);
%     C5(i,:) = log10(abs( sol(end,:)' - zf)) ;
% end
%     


figure(2)
subplot(1,2,1)
plot(log10(nfe),C1(:,1))
hold on 
plot(log10(nfe),C2(:,1))
hold on
plot(log10(nfe),C3(:,1))
hold on
plot(log10(nfe),C4(:,1))
% hold on
% plot(log10(nfe),C5(:,1))
xlabel('log10(fe)') 
ylabel('log10(erreur pour z1)') 


%figure(3)
subplot(1,2,2)
plot(log10(nfe),C1(:,2))
hold on 
plot(log10(nfe),C2(:,2))
hold on
plot(log10(nfe),C3(:,2))
hold on
plot(log10(nfe),C4(:,2))
% hold on
% plot(log10(nfe),C5(:,2))

xlabel('log10(fe)') 
ylabel('log10(erreur pour z2)') 

legend('Euler','Runge', 'Heun', 'RK4');%, 'ODE45');





% 
% Résolution avec ode45
% ---------------------
% options = odeset('Stats','on');
% [T,Z] = ode45(hvfun,[t0,tf],z0,options);
%  figure(1)
%  plot_sol(T,Z)
% legend('Euler','Runge', 'Heun', 'RK4','ode45');


   





% -----------------------------------------------------------------------
% deuxieme membre de l'edo
% ------------------------
function zpoint = hvfun_pb1_energie(t,z)
    zpoint = [-z(1)+z(2) ; z(2)];
     
end

%
% Solution exacte de l'edo
% ------------------------
function zf = exphvfun_pb1_energie(t,z0)
%
    zf = [(z0(2)/2)*exp(t)+(z0(1)-z0(2)/2)*exp(-t); z0(2)*exp(t)];
end
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
    ylabel('x_1(t)')
    subplot(2,2,2)
    hold on
    plot(T,Z(:,2))
    xlabel('t')
    ylabel('x_2(t)')
    subplot(2,2,3:4)
    hold on
    plot(Z(:,1),Z(:,2))
    xlabel('x_1(t)')
    ylabel('x_2(t)')
end




