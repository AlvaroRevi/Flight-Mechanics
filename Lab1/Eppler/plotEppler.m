clear; clc; close all; 
%% Data processing
% Headers of the raw data from XFLR5
Headers = { 'alpha','Beta','CL','CDi','CDv','CD','CY', 'Cl','Cm','Cn','Cni','QInf','XCP'};

% Legend for the different cases of study 
leyenda = {'0.026m','0.025m','0.024m','0.02m','0.03m','0.01m'};

% Read one by one each of the .txt of the directory and store the data 
F = dir('*.txt');
for ii = 1:length(F)
    data{ii} =readtable(F(ii).name,'HeaderLines',5);                
    data{ii}.Properties.VariableNames = Headers;
    
    alpha(:,ii) = data{ii}.alpha;
    Cm(:,ii) = data{ii}.Cm;
    CL(:,ii) = data{ii}.CL; 

end

%% Plots 
figure(1)
hold on 
plot(alpha,Cm)
yline(0,'--')
grid minor
axis square
xlim([-4,7.5])
xlabel('$\alpha$ [rad]','Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
legend(leyenda,'Interpreter','latex')
title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')

figure(2)
hold on 
plot(alpha,CL)
yline(0,'--')
grid minor
axis square
xlim([-4,7.5])
xlabel('$\alpha$ [rad]','Interpreter','latex')
ylabel('$C_{L}$','Interpreter','latex')
legend(leyenda,'Interpreter','latex')
title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')

figure(3)
hold on 
plot(CL,Cm)
yline(0,'--')
xlim([-0.5,0.6])
grid minor
axis square
xlabel('$C_{L}$' ,'Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
legend(leyenda,'Interpreter','latex')
title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')
