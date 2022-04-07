clear; clc; close all; 
%% Global parameters 
AR = 15.172;
AR_elev = [4.8, 2.4, 4.8, 2.4];
S_elev = [0.01,0.02,0.04,0.005];
m = [0.72, 0.735, 0.78, 0.72] ;  
g = 9.806; 
W = m*g;
S = [0.34, 0.35, 0.38, 0.34];
rho = 1.225; 
deltaH = 20; 
%% Data processing
% Headers of the raw data from XFLR5
Headers = { 'alpha','Beta','CL','CDi','CDv','CD','CY', 'Cl','Cm','Cn','Cni','QInf','XCP'};

% Legend for the different cases of study 

leyenda = {'Nominal Case','$S_{elev}=0.02m^{2}$, AR = 2.4 ','$S_{elev}=0.04m^{2}$, AR = 4.8','$S_{elev}=0.005m^{2}$, AR = 2.4'};


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

V_Nom = sqrt(2*W(1)/(rho*S(1)*Cl_zero_moment_Nom));


R_Nom = deltaH.*CL_Nom(CL_Nom>=0)./CD_Nom(CL_Nom>=0);

rate_descent_Nom = - sqrt(2*W(1)/(rho*S(1))).*CD_Nom(CL_Nom>=0)./(CL_Nom(CL_Nom>=0).^(3/2)); 

Endurance_Nom = -deltaH./rate_descent_Nom;

% Trim condition

CD_zero_moment_Nom = interp1(alpha_Nom-alpha_0_Nom,CD_Nom,alpha_zero_torque_Nom);


R_Trim_Nom = deltaH*Cl_zero_moment_Nom/CD_zero_moment_Nom;

E_Trim_Nom = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_Nom/(Cl_zero_moment_Nom^(3/2)));

%% Elevator surface = 0.2m^2 AR = 2.4

data_S2 =readtable('S2_T1-20_1 m_s-VLM2.txt','HeaderLines',5);                
data_S2.Properties.VariableNames = Headers;

alpha_S2 = data_S2.alpha;
Cm_S2 = data_S2.Cm;
CL_S2 = data_S2.CL; 
CD_S2 = data_S2.CD; 

for i=1:length(alpha_S2)
  alpha_0_S2 = interp1(CL_S2,alpha_S2,0);
end

alpha_zero_torque_S2 = interp1(Cm_S2,alpha_S2-alpha_0_S2,0);

Cl_zero_moment_S2 = interp1(alpha_S2-alpha_0_S2,CL_S2,alpha_zero_torque_S2);

V_S2 = sqrt(2*W(2)/(rho*S(2)*Cl_zero_moment_S2));


R_S2 = deltaH.*CL_S2(CL_S2>=0)./CD_S2(CL_S2>=0);

rate_descent_S2 = - sqrt(2*W(2)/(rho*S(2))).*CD_S2(CL_S2>=0)./(CL_S2(CL_S2>=0).^(3/2)); 

Endurance_S2 = -deltaH./rate_descent_S2;

CD_zero_moment_S2 = interp1(alpha_S2-alpha_0_S2,CD_S2,alpha_zero_torque_S2);

R_Trim_S2 = deltaH*Cl_zero_moment_S2/CD_zero_moment_S2;

E_Trim_S2 = deltaH/(sqrt(2*W(2)/(rho*S(2))).*CD_zero_moment_S2/(Cl_zero_moment_S2^(3/2)));

%% Surface elevator = 0.05m^2 , AR = 4.8
data_S4 =readtable('S4_T1-27_1 m_s-VLM2.txt','HeaderLines',5);                
data_S4.Properties.VariableNames = Headers;

alpha_S4 = data_S4.alpha;
Cm_S4 = data_S4.Cm;
CL_S4 = data_S4.CL; 
CD_S4 = data_S4.CD; 

for i=1:length(alpha_S4)
  alpha_0_S4 = interp1(CL_S4,alpha_S4,0);
end

alpha_zero_torque_S4 = interp1(Cm_S4,alpha_S4-alpha_0_S4,0);

Cl_zero_moment_S4 = interp1(alpha_S4-alpha_0_S4,CL_S4,alpha_zero_torque_S4);

V_S4 = sqrt(2*W(3)/(rho*S(3)*Cl_zero_moment_S4));


R_S4 = deltaH.*CL_S4(CL_S4>=0)./CD_S4(CL_S4>=0);

rate_descent_S4 = - sqrt(2*W(3)/(rho*S(3))).*CD_S4(CL_S4>=0)./(CL_S4(CL_S4>=0).^(3/2)); 

Endurance_S4 = -deltaH./rate_descent_S4;

% Trim condition 
CD_zero_moment_S4 = interp1(alpha_S4-alpha_0_S4,CD_S4,alpha_zero_torque_S4);

R_Trim_S4 = deltaH*Cl_zero_moment_S4/CD_zero_moment_S4;

E_Trim_S4 = deltaH/(sqrt(2*W(3)/(rho*S(3))).*CD_zero_moment_S4/(Cl_zero_moment_S4^(3/2)));

%% Surface elevator = 0.01m^2  AR = 2.4

data_S05 =readtable('S05_T1-19_6 m_s-VLM2.txt','HeaderLines',5);                
data_S05.Properties.VariableNames = Headers;

alpha_S05 = data_S05.alpha;
Cm_S05 = data_S05.Cm;
CL_S05 = data_S05.CL; 
CD_S05 = data_S05.CD; 

for i=1:length(alpha_S05)
  alpha_0_S05 = interp1(CL_S05,alpha_S05,0);
end

alpha_zero_torque_S05 = interp1(Cm_S05,alpha_S05-alpha_0_S05,0);

Cl_zero_moment_S05 = interp1(alpha_S05-alpha_0_S05,CL_S05,alpha_zero_torque_S05);

V_S05 = sqrt(2*W(4)/(rho*S(4)*Cl_zero_moment_S05));


R_S05 = deltaH.*CL_S05(CL_S05>=0)./CD_S05(CL_S05>=0);

rate_descent_S05 = - sqrt(2*W(4)/(rho*S(4))).*CD_S05(CL_S05>=0)./(CL_S05(CL_S05>=0).^(3/2)); 

Endurance_S05 = -deltaH./rate_descent_S05;

% Trim condition 
CD_zero_moment_S05 = interp1(alpha_S05-alpha_0_S05,CD_S05,alpha_zero_torque_S05);

R_Trim_S05 = deltaH*Cl_zero_moment_S05/CD_zero_moment_S05;

E_Trim_S05 = deltaH/(sqrt(2*W(4)/(rho*S(4))).*CD_zero_moment_S05/(Cl_zero_moment_S05^(3/2)));

%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_S2-alpha_0_S2,Cm_S2)
plot(alpha_S4-alpha_0_S4,Cm_S4)
plot(alpha_S05-alpha_0_S05,Cm_S05)
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
plot(alpha_S2-alpha_0_S2,CL_S2)
plot(alpha_S4-alpha_0_S4,CL_S4)
plot(alpha_S05-alpha_0_S05,CL_S05)
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
plot(Cm_S2,CL_S2)
plot(Cm_S4,CL_S4)
plot(Cm_S05,CL_S05)
yline(0,'--')
xlim([-0.5,0.6])
grid minor
axis square
xlabel('$C_{L}$' ,'Interpreter','latex')
ylabel('$C_{m}$','Interpreter','latex')
legend(leyenda,'Interpreter','latex')
title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')

figure(4)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,Endurance_Nom,'b--','LineWidth',1)
plot(alpha_S2(CL_S2>=0)-alpha_0_S2,Endurance_S2,'r-','LineWidth',1)
plot(alpha_S4(CL_S4>=0)-alpha_0_S4,Endurance_S4,'m-','LineWidth',1)
plot(alpha_S05(CL_S05>=0)-alpha_0_S05,Endurance_S05,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,E_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_S2,E_Trim_S2,'r.','MarkerSize',20)
plot(alpha_zero_torque_S4,E_Trim_S4,'m.','MarkerSize',20)
plot(alpha_zero_torque_S05,E_Trim_S05,'k.','MarkerSize',20)

yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex','location','best')
title('Aircraft performance for different elevator size','Interpreter','latex','FontSize',14)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_S2(CL_S2>=0)-alpha_0_S2,R_S2,'r-','LineWidth',1)
plot(alpha_S4(CL_S4>=0)-alpha_0_S4,R_S4,'m-','LineWidth',1)
plot(alpha_S05(CL_S05>=0)-alpha_0_S05,R_S05,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_S2,R_Trim_S2,'r.','MarkerSize',20)
plot(alpha_zero_torque_S4,R_Trim_S4,'m.','MarkerSize',20)
plot(alpha_zero_torque_S05,R_Trim_S05,'k.','MarkerSize',20)

yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Range [m]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex','Location','best')
title('Aircraft performance for different elevator size','Interpreter','latex','FontSize',14)

figure(6)
hold on
grid minor
axis square
plot3(AR_elev,S_elev, [V_Nom,V_S2,V_S4,V_S05],'r.','MarkerSize',20)
xlabel('Aspect Ratio','Interpreter','latex','FontSize',14)
ylabel('Elevator Surface [m$^{2}$]','Interpreter','latex','FontSize',14)
zlabel('Trim Velocity [m/s]','Interpreter','latex','FontSize',14)
title('Trim condition for different elevator sizes','Interpreter','latex','FontSize',18)
