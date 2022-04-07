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

leyenda = {'Nominal Case','$U_{\infty}$=30m/s','$U_{\infty}$=40m/s','$U_{\infty}$=16m/s'};


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

E_Trim_Nom = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_Nom/(Cl_zero_moment_Nom^(3/2)));


%% U_inf = 30 m/s 
data_30 =readtable('T1-30_0 m_s-VLM2.txt','HeaderLines',5);                
data_30.Properties.VariableNames = Headers;

alpha_30 = data_30.alpha;
Cm_30 = data_30.Cm;
CL_30 = data_30.CL; 
CD_30 = data_30.CD; 

for i=1:length(alpha_30)
  alpha_0_30 = interp1(CL_30,alpha_30,0);
end

alpha_zero_torque_30 = interp1(Cm_30,alpha_30-alpha_0_30,0);

Cl_zero_moment_30 = interp1(alpha_30-alpha_0_30,CL_30,alpha_zero_torque_30);

V_30 = sqrt(2*W/(rho*S*Cl_zero_moment_30));


R_30 = deltaH.*CL_30(CL_30>=0)./CD_30(CL_30>=0);

rate_descent_30 = - sqrt(2*W/(rho*S)).*CD_30(CL_30>=0)./(CL_30(CL_30>=0).^(3/2)); 

Endurance_30 = -deltaH./rate_descent_30;

% Trim condition

CD_zero_moment_30 = interp1(alpha_30-alpha_0_30,CD_30,alpha_zero_torque_30);

R_Trim_30 = deltaH*Cl_zero_moment_30/CD_zero_moment_30;

E_Trim_30 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_30/(Cl_zero_moment_30^(3/2)));

%% U_inf = 17

data_40 =readtable('T1-40_0 m_s-VLM2.txt','HeaderLines',5);                
data_40.Properties.VariableNames = Headers;

alpha_40 = data_40.alpha;
Cm_40 = data_40.Cm;
CL_40 = data_40.CL; 
CD_40 = data_40.CD; 

for i=1:length(alpha_40)
  alpha_0_40 = interp1(CL_40,alpha_40,0);
end

alpha_zero_torque_40 = interp1(Cm_40,alpha_40-alpha_0_40,0);

Cl_zero_moment_40 = interp1(alpha_40-alpha_0_40,CL_40,alpha_zero_torque_40);

V_40 = sqrt(2*W/(rho*S*Cl_zero_moment_40));

R_40 = deltaH.*CL_40(CL_40>=0)./CD_40(CL_40>=0);

rate_descent_40 = - sqrt(2*W/(rho*S)).*CD_40(CL_40>=0)./(CL_40(CL_40>=0).^(3/2)); 

Endurance_40 = -deltaH./rate_descent_40;

% Trim condition

CD_zero_moment_40 = interp1(alpha_40-alpha_0_40,CD_40,alpha_zero_torque_40);

R_Trim_40 = deltaH*Cl_zero_moment_40/CD_zero_moment_40;

E_Trim_40 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_40/(Cl_zero_moment_40^(3/2)));

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

% Trim condition

CD_zero_moment_16 = interp1(alpha_16-alpha_0_16,CD_16,alpha_zero_torque_16);

R_Trim_16 = deltaH*Cl_zero_moment_16/CD_zero_moment_16;

E_Trim_16 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_16/(Cl_zero_moment_16^(3/2)));

%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_30-alpha_0_30,Cm_30)
plot(alpha_40-alpha_0_40,Cm_40)
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
plot(alpha_30-alpha_0_30,CL_30)
plot(alpha_40-alpha_0_40,CL_40)
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
plot(Cm_30,CL_30)
plot(Cm_40,CL_40)
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
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,Endurance_Nom,'b--','LineWidth',1)
plot(alpha_40(CL_40>=0)-alpha_0_40,Endurance_40,'r-','LineWidth',1)
plot(alpha_30(CL_30>=0)-alpha_0_30,Endurance_30,'m-','LineWidth',1)
plot(alpha_16(CL_16>=0)-alpha_0_16,Endurance_16,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,E_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_40,E_Trim_40,'r.','MarkerSize',20)
plot(alpha_zero_torque_30,E_Trim_30,'m.','MarkerSize',20)
plot(alpha_zero_torque_16,E_Trim_16,'k.','MarkerSize',20)
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
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_40(CL_40>=0)-alpha_0_40,R_40,'r-','LineWidth',1)
plot(alpha_30(CL_30>=0)-alpha_0_30,R_30,'m-','LineWidth',1)
plot(alpha_16(CL_16>=0)-alpha_0_16,R_16,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_40,R_Trim_40,'r.','MarkerSize',20)
plot(alpha_zero_torque_30,R_Trim_30,'m.','MarkerSize',20)
plot(alpha_zero_torque_16,R_Trim_16,'k.','MarkerSize',20)
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
% plot([-6,-4,-1,1], [V_30,V_40,V_Nom,V_16],'r.-','MarkerSize',10,'LineWidth',2)
% xlabel('Twist at the tip [$^{\circ}$]','Interpreter','latex','FontSize',18)
% ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
% title('Trim condition for different twist angles at the tip','Interpreter','latex','FontSize',18)
