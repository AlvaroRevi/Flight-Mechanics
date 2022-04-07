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

leyenda = {'Nominal Case','$X_{elev}= 0.5m$','$X_{elev}= 0.3m$','$X_{elev}= -0.1m$'};


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

% Trim condition

CD_zero_moment_Nom = interp1(alpha_Nom-alpha_0_Nom,CD_Nom,alpha_zero_torque_Nom);

R_Trim_Nom = deltaH*Cl_zero_moment_Nom/CD_zero_moment_Nom;

E_Trim_Nom = deltaH/(sqrt(2*W/(rho*S)).*CD_zero_moment_Nom/(Cl_zero_moment_Nom^(3/2)));

%% Elevator distance = 0.5m

data_5 =readtable('ED_5_T1-19_3 m_s-VLM2.txt','HeaderLines',5);                
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

V_5 = sqrt(2*W/(rho*S*Cl_zero_moment_5));


R_5 = deltaH.*CL_5(CL_5>=0)./CD_5(CL_5>=0);

rate_descent_5 = - sqrt(2*W/(rho*S)).*CD_5(CL_5>=0)./(CL_5(CL_5>=0).^(3/2)); 

Endurance_5 = -deltaH./rate_descent_5;

CD_zero_moment_5 = interp1(alpha_5-alpha_0_5,CD_5,alpha_zero_torque_5);

R_Trim_5 = deltaH*Cl_zero_moment_5/CD_zero_moment_5;

E_Trim_5 = deltaH/(sqrt(2*W/(rho*S)).*CD_zero_moment_5/(Cl_zero_moment_5^(3/2)));

%% Elevator distance = 0.3 

data_3 =readtable('ED_3_T1-19_24 m_s-VLM2.txt','HeaderLines',5);                
data_3.Properties.VariableNames = Headers;

alpha_3 = data_3.alpha;
Cm_3 = data_3.Cm;
CL_3 = data_3.CL; 
CD_3 = data_3.CD; 

for i=1:length(alpha_3)
  alpha_0_3 = interp1(CL_3,alpha_3,0,'pchip');
end

alpha_zero_torque_3 = interp1(Cm_3,alpha_3-alpha_0_3,0,'pchip');

Cl_zero_moment_3 = interp1(alpha_3-alpha_0_3,CL_3,alpha_zero_torque_3,'pchip');

V_3 = sqrt(2*W/(rho*S*Cl_zero_moment_3));


R_3 = deltaH.*CL_3(CL_3>=0)./CD_3(CL_3>=0);

rate_descent_3 = - sqrt(2*W/(rho*S)).*CD_3(CL_3>=0)./(CL_3(CL_3>=0).^(3/2)); 

Endurance_3 = -deltaH./rate_descent_3;

CD_zero_moment_3 = interp1(alpha_3-alpha_0_3,CD_3,alpha_zero_torque_3);

R_Trim_3 = deltaH*Cl_zero_moment_3/CD_zero_moment_3;

E_Trim_3 = deltaH/(sqrt(2*W/(rho*S)).*CD_zero_moment_3/(Cl_zero_moment_3^(3/2)));

%% Elevator distance = -0.1m

data_1 =readtable('ED_m1_T1-20_85 m_s-VLM2.txt','HeaderLines',5);                
data_1.Properties.VariableNames = Headers;

alpha_1 = data_1.alpha;
Cm_1 = data_1.Cm;
CL_1 = data_1.CL; 
CD_1 = data_1.CD; 

for i=1:length(alpha_1)
  alpha_0_1 = interp1(CL_1,alpha_1,0);
end

alpha_zero_torque_1 = interp1(Cm_1,alpha_1-alpha_0_1,0);

Cl_zero_moment_1 = interp1(alpha_1-alpha_0_1,CL_1,alpha_zero_torque_1);

V_1 = sqrt(2*W/(rho*S*Cl_zero_moment_1));


R_1 = deltaH.*CL_1(CL_1>=0)./CD_1(CL_1>=0);

rate_descent_1 = - sqrt(2*W/(rho*S)).*CD_1(CL_1>=0)./(CL_1(CL_1>=0).^(3/2)); 

Endurance_1 = -deltaH./rate_descent_1;

CD_zero_moment_1 = interp1(alpha_1-alpha_0_1,CD_1,alpha_zero_torque_1);

R_Trim_1 = deltaH*Cl_zero_moment_1/CD_zero_moment_1;

E_Trim_1 = deltaH/(sqrt(2*W/(rho*S)).*CD_zero_moment_1/(Cl_zero_moment_1^(3/2)));

%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_5-alpha_0_5,Cm_5)
plot(alpha_3-alpha_0_3,Cm_3)
plot(alpha_1-alpha_0_1,Cm_1)
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
plot(alpha_5-alpha_0_5,CL_5)
plot(alpha_3-alpha_0_3,CL_3)
plot(alpha_1-alpha_0_1,CL_1)
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
plot(Cm_5,CL_5)
plot(Cm_3,CL_3)
plot(Cm_1,CL_1)
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
plot(alpha_5(CL_5>=0)-alpha_0_5,Endurance_5,'r-','LineWidth',1)
plot(alpha_3(CL_3>=0)-alpha_0_3,Endurance_3,'m-','LineWidth',1)
plot(alpha_1(CL_1>=0)-alpha_0_1,Endurance_1,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,E_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_5,E_Trim_5,'r.','MarkerSize',20)
plot(alpha_zero_torque_3,E_Trim_3,'m.','MarkerSize',20)
plot(alpha_zero_torque_1,E_Trim_1,'k.','MarkerSize',20)

yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [rad]','Interpreter','latex','FontSize',14)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for elevator position','Interpreter','latex','FontSize',14)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_5(CL_5>=0)-alpha_0_5,R_5,'r-','LineWidth',1)
plot(alpha_3(CL_3>=0)-alpha_0_3,R_3,'m-','LineWidth',1)
plot(alpha_1(CL_1>=0)-alpha_0_1,R_1,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_5,R_Trim_5,'r.','MarkerSize',20)
plot(alpha_zero_torque_3,R_Trim_3,'m.','MarkerSize',20)
plot(alpha_zero_torque_1,R_Trim_1,'k.','MarkerSize',20)

yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Range [m]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different elevator position','Interpreter','latex','FontSize',14)

figure(6)
hold on
grid minor
axis square
plot([-0.1,0.3,0.5,0.6], [V_1,V_3,V_5,V_Nom],'r.-','MarkerSize',10,'LineWidth',2)
xlabel('Elevator position in body basis [m]','Interpreter','latex','FontSize',18)
ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
title('Trim condition for different elevator position','Interpreter','latex','FontSize',18)
