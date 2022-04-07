clear; clc; close all; 
%% Global parameters 
AR = [15.172, 6.89, 11.03, 27.59];
m = [0.72, 0.45, 0.675, 1.11];  
g = 9.806; 
W = m*g;
S = [0.34, 0.145,0.29, 0.58];
rho = 1.225; 
deltaH = 20; 
%% Data processing
% Headers of the raw data from XFLR5
Headers = { 'alpha','Beta','CL','CDi','CDv','CD','CY', 'Cl','Cm','Cn','Cni','QInf','XCP'};

% Legend for the different cases of study 

leyenda = {'Nominal Case','AR = 6.89','AR = 11.03','AR = 27.59'};


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

%% AR = 6.89 
data_6 =readtable('AR_6_89_T1-39_01_m_s-VLM2.txt','HeaderLines',5);                
data_6.Properties.VariableNames = Headers;

alpha_6 = data_6.alpha;
Cm_6 = data_6.Cm;
CL_6 = data_6.CL; 
CD_6 = data_6.CD; 

for i=1:length(alpha_6)
  alpha_0_6 = interp1(CL_6,alpha_6,0);
end

alpha_zero_torque_6 = interp1(Cm_6,alpha_6-alpha_0_6,0);

Cl_zero_moment_6 = interp1(alpha_6-alpha_0_6,CL_6,alpha_zero_torque_6);

V_6 = sqrt(2*W(2)/(rho*S(2)*Cl_zero_moment_6));

R_6 = deltaH.*CL_6(CL_6>=0)./CD_6(CL_6>=0);

rate_descent_6 = - sqrt(2*W(2)/(rho*S(2))).*CD_6(CL_6>=0)./(CL_6(CL_6>=0).^(3/2)); 

Endurance_6 = -deltaH./rate_descent_6;


%% AR = 11.03
data_11 =readtable('AR_11_03_T1-26_44 m_s-VLM2.txt','HeaderLines',5);                
data_11.Properties.VariableNames = Headers;

alpha_11 = data_11.alpha;
Cm_11 = data_11.Cm;
CL_11 = data_11.CL; 
CD_11 = data_11.CD; 

for i=1:length(alpha_11)
  alpha_0_11 = interp1(CL_11,alpha_11,0);
end

alpha_zero_torque_11 = interp1(Cm_11,alpha_11-alpha_0_11,0);

Cl_zero_moment_11 = interp1(alpha_11-alpha_0_11,CL_11,alpha_zero_torque_11);

V_11 = sqrt(2*W(3)/(rho*S(3)*Cl_zero_moment_11));

R_11 = deltaH.*CL_11(CL_11>=0)./CD_11(CL_11>=0);

rate_descent_11 = - sqrt(2*W(3)/(rho*S(3))).*CD_11(CL_11>=0)./(CL_11(CL_11>=0).^(3/2)); 

Endurance_11 = -deltaH./rate_descent_11;

%% AR = 27.59
data_27 =readtable('AR_27_586_T1-15 m_s-VLM2.txt','HeaderLines',5);                
data_27.Properties.VariableNames = Headers;

alpha_27 = data_27.alpha;
Cm_27 = data_27.Cm;
CL_27 = data_27.CL; 
CD_27 = data_27.CD; 

for i=1:length(alpha_27)
  alpha_0_27 = interp1(CL_27,alpha_27,0);
end

alpha_zero_torque_27 = interp1(Cm_27,alpha_27-alpha_0_27,0);

Cl_zero_moment_27 = interp1(alpha_27-alpha_0_27,CL_27,alpha_zero_torque_27);

V_27 = sqrt(2*W(4)/(rho*S(4)*Cl_zero_moment_27));

%%%%%%% NOTE: The trim velocity is found at approx 13.22 m/s. However, the
%%%%%%% minimum velocity allowed by XFLR5 is 15.5 m/s in order to avoid 
%%%%%%% convergence errors. Hence, it is assumed 13.22 as the trim velocity
%%%%%%% in this case of study, as convergence has not been reached. 

R_27 = deltaH.*CL_27(CL_27>=0)./CD_27(CL_27>=0);

rate_descent_27 = - sqrt(2*W(1)/(rho*S(1))).*CD_27(CL_27>=0)./(CL_27(CL_27>=0).^(3/2)); 

Endurance_27 = -deltaH./rate_descent_27;


%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_6-alpha_0_6,Cm_6)
plot(alpha_11-alpha_0_11,Cm_11)
plot(alpha_27-alpha_0_27,Cm_27)
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
plot(alpha_6-alpha_0_6,CL_6)
plot(alpha_11-alpha_0_11,CL_11)
plot(alpha_27-alpha_0_27,CL_27)
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
plot(Cm_6,CL_6)
plot(Cm_11,CL_11)
plot(Cm_27,CL_27)
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
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,Endurance_Nom,'--','LineWidth',1)
plot(alpha_6(CL_6>=0)-alpha_0_6,Endurance_6,'LineWidth',1)
plot(alpha_11(CL_11>=0)-alpha_0_11,Endurance_11,'LineWidth',1)
plot(alpha_27(CL_27>=0)-alpha_0_27,Endurance_27,'LineWidth',1)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',18)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',18)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different aspect ratio','Interpreter','latex','FontSize',18)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_6(CL_6>=0)-alpha_0_6,R_6,'r-','LineWidth',1)
plot(alpha_11(CL_11>=0)-alpha_0_11,R_11,'g-','LineWidth',1)
plot(alpha_27(CL_27>=0)-alpha_0_27,R_27,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b*','MarkerSize',10)

yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',18)
ylabel('$Range [m]$','Interpreter','latex','FontSize',18)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different aspect ratio','Interpreter','latex','FontSize',18)

figure(6)
hold on
grid minor
axis square
plot([6.89, 11.03, 15.172 27.59], [V_6,V_11,V_Nom,V_27],'r.-','MarkerSize',10,'LineWidth',2)
xlabel('Aspect Ratio','Interpreter','latex','FontSize',18)
ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
title('Trim condition for different aspect ratio','Interpreter','latex','FontSize',18)
