clear; clc; close all; 
%% Global parameters 
AR = 15.172;
m = 0.72;  
g = 9.806; 
W = m*g;
S = 0.34;
rho = 1.225; 
deltaH = 20; 
%% Data processing
% Headers of the raw data from XFLR5
Headers = { 'alpha','Beta','CL','CDi','CDv','CD','CY', 'Cl','Cm','Cn','Cni','QInf','XCP'};

% Legend for the different cases of study 

leyenda = {'Nominal Case','$U_{\infty}$=18m/s','$U_{\infty}$=17m/s','$U_{\infty}$=16m/s'};


%% Nominal case 

data_Nom =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_Nom.Properties.VariableNames = Headers;

alpha_Nom = data_Nom.alpha;
Cm_Nom = data_Nom.Cm;
CL_Nom = data_Nom.CL; 
CD_Nom = data_Nom.CD; 

for i=1:length(alpha_Nom)
  alpha_0_Nom = interp1(CL_Nom,alpha_Nom,0);
end

alpha_zero_torque_Nom = interp1(Cm_Nom,alpha_Nom-alpha_0_Nom,0);

Cl_zero_moment_Nom = interp1(alpha_Nom-alpha_0_Nom,CL_Nom,alpha_zero_torque_Nom);

V_Nom = sqrt(2*W/(rho*S*Cl_zero_moment_Nom));


R_Nom = deltaH.*CL_Nom(CL_Nom>=0)./CD_Nom(CL_Nom>=0);

rate_descent_Nom = - sqrt(2*W/(rho*S)).*CD_Nom(CL_Nom>=0)./(CL_Nom(CL_Nom>=0).^(3/2)); 

Endurance_Nom = -deltaH./rate_descent_Nom;

%% U_inf = 18 m/s 
data_18 =readtable('T1-18_0 m_s-VLM2.txt','HeaderLines',5);                
data_18.Properties.VariableNames = Headers;

alpha_18 = data_18.alpha;
Cm_18 = data_18.Cm;
CL_18 = data_18.CL; 
CD_18 = data_18.CD; 

for i=1:length(alpha_18)
  alpha_0_18 = interp1(CL_18,alpha_18,0);
end

alpha_zero_torque_18 = interp1(Cm_18,alpha_18-alpha_0_18,0);

Cl_zero_moment_18 = interp1(alpha_18-alpha_0_18,CL_18,alpha_zero_torque_18);

V_18 = sqrt(2*W/(rho*S*Cl_zero_moment_18));


R_18 = deltaH.*CL_18(CL_18>=0)./CD_18(CL_18>=0);

rate_descent_18 = - sqrt(2*W/(rho*S)).*CD_18(CL_18>=0)./(CL_18(CL_18>=0).^(3/2)); 

Endurance_18 = -deltaH./rate_descent_18;


%% U_inf = 17

data_17 =readtable('T1-17_0 m_s-VLM2.txt','HeaderLines',5);                
data_17.Properties.VariableNames = Headers;

alpha_17 = data_17.alpha;
Cm_17 = data_17.Cm;
CL_17 = data_17.CL; 
CD_17 = data_17.CD; 

for i=1:length(alpha_17)
  alpha_0_17 = interp1(CL_17,alpha_17,0);
end

alpha_zero_torque_17 = interp1(Cm_17,alpha_17-alpha_0_17,0);

Cl_zero_moment_17 = interp1(alpha_17-alpha_0_17,CL_17,alpha_zero_torque_17);

V_17 = sqrt(2*W/(rho*S*Cl_zero_moment_17));


R_17 = deltaH.*CL_17(CL_17>=0)./CD_17(CL_17>=0);

rate_descent_17 = - sqrt(2*W/(rho*S)).*CD_17(CL_17>=0)./(CL_17(CL_17>=0).^(3/2)); 

Endurance_17 = -deltaH./rate_descent_17;

%% U_inf = 16 m/s 

data_16 =readtable('T1-16_0 m_s-VLM2.txt','HeaderLines',5);                
data_16.Properties.VariableNames = Headers;

alpha_16 = data_16.alpha;
Cm_16 = data_16.Cm;
CL_16 = data_16.CL; 
CD_16 = data_16.CD; 

for i=1:length(alpha_16)
  alpha_0_16 = interp1(CL_16,alpha_16,0);
end

alpha_zero_torque_16 = interp1(Cm_16,alpha_16-alpha_0_16,0);

Cl_zero_moment_16 = interp1(alpha_16-alpha_0_16,CL_16,alpha_zero_torque_16);

V_16 = sqrt(2*W/(rho*S*Cl_zero_moment_16));


R_16 = deltaH.*CL_16(CL_16>=0)./CD_16(CL_16>=0);

rate_descent_16 = - sqrt(2*W/(rho*S)).*CD_16(CL_16>=0)./(CL_16(CL_16>=0).^(3/2)); 

Endurance_16 = -deltaH./rate_descent_16;

%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_18-alpha_0_18,Cm_18)
plot(alpha_17-alpha_0_17,Cm_17)
plot(alpha_16-alpha_0_16,Cm_16)
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
plot(alpha_Nom-alpha_0_Nom,CL_Nom)
plot(alpha_18-alpha_0_18,CL_18)
plot(alpha_17-alpha_0_17,CL_17)
plot(alpha_16-alpha_0_16,CL_16)
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
plot(Cm_Nom,CL_Nom)
plot(Cm_18,CL_18)
plot(Cm_17,CL_17)
plot(Cm_16,CL_16)
yline(0,'--')
xlim([-0.5,0.6])
grid minor
axis square
xlabel('$C_{L}$' ,'Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
legend(leyenda,'Interpreter','latex')
title('Analysis of the aircraft performance for different headwinds','Interpreter','latex')

figure(4)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,Endurance_Nom,'--','LineWidth',1)
plot(alpha_18(CL_18>=0)-alpha_0_18,Endurance_18,'LineWidth',1)
plot(alpha_17(CL_17>=0)-alpha_0_17,Endurance_17,'LineWidth',1)
plot(alpha_16(CL_16>=0)-alpha_0_16,Endurance_16,'LineWidth',1)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',18)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',18)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different headwinds','Interpreter','latex','FontSize',18)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'--','LineWidth',1)
plot(alpha_18(CL_18>=0)-alpha_0_18,R_18,'LineWidth',1)
plot(alpha_17(CL_17>=0)-alpha_0_17,R_17,'LineWidth',1)
plot(alpha_16(CL_16>=0)-alpha_0_16,R_16,'LineWidth',1)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',18)
ylabel('$Range [m]$','Interpreter','latex','FontSize',18)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different headwinds','Interpreter','latex','FontSize',18)

% figure(6)
% hold on
% grid minor
% axis square
% plot([-6,-4,-1,1], [V_18,V_17,V_Nom,V_16],'r.-','MarkerSize',10,'LineWidth',2)
% xlabel('Twist at the tip [$^{\circ}$]','Interpreter','latex','FontSize',18)
% ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
% title('Trim condition for different twist angles at the tip','Interpreter','latex','FontSize',18)
