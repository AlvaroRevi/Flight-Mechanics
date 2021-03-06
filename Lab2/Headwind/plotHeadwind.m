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

leyenda = {'Nominal Case','$V_{g}$=10m/s','$V_{g}$=15m/s','$V_{g}$=5m/s'};


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

% Endurance_Nom = -deltaH./rate_descent_Nom;

Endurance_Nom = R_Nom/V_Nom;

% Trim condition

CD_zero_moment_Nom = interp1(alpha_Nom-alpha_0_Nom,CD_Nom,alpha_zero_torque_Nom);

R_Trim_Nom = deltaH*Cl_zero_moment_Nom/CD_zero_moment_Nom;

%E_Trim_Nom = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_Nom/(Cl_zero_moment_Nom^(3/2)));
E_Trim_Nom = R_Trim_Nom/V_Nom;  

%% U_inf = 15 m/s 
data_15 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_15.Properties.VariableNames = Headers;

alpha_15 = data_15.alpha;
Cm_15 = data_15.Cm;
CL_15 = data_15.CL; 
CD_15 = data_15.CD; 

for i=1:length(alpha_15)
  alpha_0_15 = interp1(CL_15,alpha_15,0);
end

alpha_zero_torque_15 = interp1(Cm_15,alpha_15-alpha_0_15,0);

Cl_zero_moment_15 = interp1(alpha_15-alpha_0_15,CL_15,alpha_zero_torque_15);

V_15 = sqrt(2*W/(rho*S*Cl_zero_moment_15)) - 5;



R_15 = deltaH.*CL_15(CL_15>=0)./CD_15(CL_15>=0);

rate_descent_15 = - sqrt(2*W/(rho*S)).*CD_15(CL_15>=0)./(CL_15(CL_15>=0).^(3/2)); 

%Endurance_15 = -deltaH./rate_descent_15;
Endurance_15 = R_15/V_15; 


% Trim condition

CD_zero_moment_15 = interp1(alpha_15-alpha_0_15,CD_15,alpha_zero_torque_15);

R_Trim_15 = deltaH*Cl_zero_moment_15/CD_zero_moment_15;

%E_Trim_15 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_15/(Cl_zero_moment_15^(3/2)));
 
E_Trim_15 = R_Trim_15/V_15;

%% U_inf = 10

data_10 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_10.Properties.VariableNames = Headers;

alpha_10 = data_10.alpha;
Cm_10 = data_10.Cm;
CL_10 = data_10.CL; 
CD_10 = data_10.CD; 

for i=1:length(alpha_10)
  alpha_0_10 = interp1(CL_10,alpha_10,0);
end

alpha_zero_torque_10 = interp1(Cm_10,alpha_10-alpha_0_10,0);

Cl_zero_moment_10 = interp1(alpha_10-alpha_0_10,CL_10,alpha_zero_torque_10);

V_10 = sqrt(2*W/(rho*S*Cl_zero_moment_10)) -10;

R_10 = deltaH.*CL_10(CL_10>=0)./CD_10(CL_10>=0);

rate_descent_10 = - sqrt(2*W/(rho*S)).*CD_10(CL_10>=0)./(CL_10(CL_10>=0).^(3/2)); 

%Endurance_10 = -deltaH./rate_descent_10;
Endurance_10 = R_10/V_10;

% Trim condition

CD_zero_moment_10 = interp1(alpha_10-alpha_0_10,CD_10,alpha_zero_torque_10);

R_Trim_10 = deltaH*Cl_zero_moment_10/CD_zero_moment_10;

%E_Trim_10 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_10/(Cl_zero_moment_10^(3/2)));
E_Trim_10 = R_Trim_10/V_10; 

%% U_inf = 5 m/s 

data_5 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_5.Properties.VariableNames = Headers;

alpha_5 = data_5.alpha;
Cm_5 = data_5.Cm;
CL_5 = data_5.CL; 
CD_5 = data_5.CD; 

for i=1:length(alpha_5)
  alpha_0_5 = interp1(CL_5,alpha_5,0);
end

alpha_zero_torque_5 = interp1(Cm_5,alpha_5-alpha_0_5,0);

Cl_zero_moment_5 = interp1(alpha_5-alpha_0_5,CL_5,alpha_zero_torque_5);

V_5 = sqrt(2*W/(rho*S*Cl_zero_moment_5)) -15;


R_5 = deltaH.*CL_5(CL_5>=0)./CD_5(CL_5>=0);

rate_descent_5 = - sqrt(2*W/(rho*S)).*CD_5(CL_5>=0)./(CL_5(CL_5>=0).^(3/2)); 

%Endurance_5 = -deltaH./rate_descent_5;
Endurance_5 = R_5/V_5; 

% Trim condition

CD_zero_moment_5 = interp1(alpha_5-alpha_0_5,CD_5,alpha_zero_torque_5);

R_Trim_5 = deltaH*Cl_zero_moment_5/CD_zero_moment_5;

%E_Tri
% m_5 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_5/(Cl_zero_moment_5^(3/2)));
E_Trim_5 = R_Trim_5/V_5; 


%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_15-alpha_0_15,Cm_15)
plot(alpha_10-alpha_0_10,Cm_10)
plot(alpha_5-alpha_0_5,Cm_5)
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
plot(alpha_15-alpha_0_15,CL_15)
plot(alpha_10-alpha_0_10,CL_10)
plot(alpha_5-alpha_0_5,CL_5)
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
plot(Cm_15,CL_15)
plot(Cm_10,CL_10)
plot(Cm_5,CL_5)
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
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,Endurance_Nom,'b--','LineWidth',1)
plot(alpha_10(CL_10>=0)-alpha_0_10,Endurance_10,'r-','LineWidth',1)
plot(alpha_15(CL_15>=0)-alpha_0_15,Endurance_15,'m-','LineWidth',1)
plot(alpha_5(CL_5>=0)-alpha_0_5,Endurance_5,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,E_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_10,E_Trim_10,'r.','MarkerSize',20)
plot(alpha_zero_torque_15,E_Trim_15,'m.','MarkerSize',20)
plot(alpha_zero_torque_5,E_Trim_5,'k.','MarkerSize',20)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different headwinds','Interpreter','latex','FontSize',14)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_10(CL_10>=0)-alpha_0_10,R_10,'r-','LineWidth',1)
plot(alpha_15(CL_15>=0)-alpha_0_15,R_15,'m-','LineWidth',1)
plot(alpha_5(CL_5>=0)-alpha_0_5,R_5,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_10,R_Trim_10,'r.','MarkerSize',20)
plot(alpha_zero_torque_15,R_Trim_15,'m.','MarkerSize',20)
plot(alpha_zero_torque_5,R_Trim_5,'k.','MarkerSize',20)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Range [m]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different headwinds','Interpreter','latex','FontSize',14)

% figure(6)
% hold on
% grid minor
% axis square
% plot([-6,-4,-1,1], [V_15,V_10,V_Nom,V_5],'r.-','MarkerSize',10,'LineWidth',2)
% xlabel('Twist at the tip [$^{\circ}$]','Interpreter','latex','FontSize',18)
% ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
% title('Trim condition for different twist angles at the tip','Interpreter','latex','FontSize',18)
