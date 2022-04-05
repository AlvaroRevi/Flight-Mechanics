clear; clc; close all; 
%% Global parameters
m = 0.72;  
g = 9.806; 
W = m*g;
S = 0.34;
rho = 1.225; 
deltaH = 20; 

%% Data processing
% Headers of the raw data from XFLR5
Headers = {'alpha','Beta','CL','CDi','CDv','CD','CY', 'Cl','Cm','Cn','Cni','QInf','XCP'};



% Read one by one each of the .txt of the directory and store the data 
F = dir('*.txt');
for ii = 1:length(F)
    data{ii} =readtable(F(ii).name,'HeaderLines',5);                
    data{ii}.Properties.VariableNames = Headers;
    
    alpha(:,ii) = data{ii}.alpha;
    Cm(:,ii) = data{ii}.Cm;
    CL(:,ii) = data{ii}.CL; 
    CD(:,ii) = data{ii}.CD;
end

%% Obtain alpha effective
C = polyfit(CL(:,1),alpha(:,1),1);
alpha_0(1) = C(2);
for i=2:length(alpha(1,:))
    alpha_0(i) = interp1(CL(:,i),alpha(:,i),0);
end

alpha_true = alpha - alpha_0; 

%% Look for the alpha that gives Cm = 0

alpha_zero_torque = interp1(Cm,alpha_true,0);

Cl_zero_moment = interp1(alpha_true,CL,alpha_zero_torque);

V = sqrt(2*W/(rho*S*Cl_zero_moment));

%% Gliding Performance 

Eff = CL./CD;

Range = deltaH.*Eff(CL>=0);

rate_descent = - sqrt(2*W/(rho*S)).*CD(CL>=0)./(CL(CL>=0).^(3/2)); 

Endurance = -deltaH./rate_descent;

%% Plots 
figure(1)
hold on 
grid minor
plot(alpha-alpha_0,Cm,'LineWidth',1)
yline(0,'--')
axis square
xlim([0,7.5])
xlabel('$\alpha$ [rad]','Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
title('Analysis of the aircraft performance','Interpreter','latex')

figure(2)
hold on 
grid minor 
plot(alpha-alpha_0,CL,'LineWidth',1)
yline(0,'--')
axis square
xlim([-3,7.5])
xlabel('$\alpha$ [rad]','Interpreter','latex')
ylabel('$C_{L}$','Interpreter','latex')
title('Analysis of the aircraft performance','Interpreter','latex')

figure(3)
hold on 
grid minor 
plot(CL,Cm,'LineWidth',1)
yline(0,'--')
xlim([-0.2,0.6])
axis square
xlabel('$C_{L}$' ,'Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
title('Analysis of the aircraft performance','Interpreter','latex')

figure(4)
hold on 
grid minor 
plot(alpha-alpha_0,Eff,'LineWidth',1)
yline(0,'--')
xlim([-4,6])
axis square
xlabel('$\alpha$' ,'Interpreter','latex')
ylabel('$C_{L}/C_{D}$','Interpreter','latex')
title('Analysis of the aircraft performance','Interpreter','latex')

figure(5)
grid minor
xlabel('Angle of attack $\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',18)
title('Glider performance','Interpreter','latex','FontSize',18)

yyaxis left 
plot(alpha_true(CL>=0),Range)
ylabel('Range [m]','Interpreter','latex','FontSize',18)

yyaxis right 
plot( alpha_true(CL>=0),Endurance)
ylabel('Endurance [s]','Interpreter','latex','FontSize',18)