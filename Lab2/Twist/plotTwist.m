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

leyenda = {'Nominal Case','$\theta$= -6$^{\circ}$','$\theta$ =-4$^{\circ}$','$\theta$ = +1$^{\circ}$'};


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

%% Twist = -6º
data_m6 =readtable('TW-6_T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_m6.Properties.VariableNames = Headers;

alpha_m6 = data_m6.alpha;
Cm_m6 = data_m6.Cm;
CL_m6 = data_m6.CL; 
CD_m6 = data_m6.CD; 

for i=1:length(alpha_m6)
  alpha_0_m6 = interp1(CL_m6,alpha_m6,0);
end

alpha_zero_torque_m6 = interp1(Cm_m6,alpha_m6-alpha_0_m6,0);

Cl_zero_moment_m6 = interp1(alpha_m6-alpha_0_m6,CL_m6,alpha_zero_torque_m6);

V_m6 = sqrt(2*W/(rho*S*Cl_zero_moment_m6));

%%%%%%% NOTE: The trim velocity is found at approx 11.4 m/s. However, the
%%%%%%% minimum velocity allowed by XFLR5 is 17 m/s in order to avoid 
%%%%%%% convergence errors. Hence, it is assumed 11.4 m/s as the trim velocity
%%%%%%% in this case of study, as convergence has not been reached. 

R_m6 = deltaH.*CL_m6(CL_m6>=0)./CD_m6(CL_m6>=0);

rate_descent_m6 = - sqrt(2*W/(rho*S)).*CD_m6(CL_m6>=0)./(CL_m6(CL_m6>=0).^(3/2)); 

Endurance_m6 = -deltaH./rate_descent_m6;

% Trim condition

CD_zero_moment_m6 = interp1(alpha_m6-alpha_0_m6,CD_m6,alpha_zero_torque_m6);

R_Trim_m6 = deltaH*Cl_zero_moment_m6/CD_zero_moment_m6;

E_Trim_m6 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_m6/(Cl_zero_moment_m6^(3/2)));
%% Twist -4º
data_m4 =readtable('TW-4_T1-17_0 m_s-VLM2.txt','HeaderLines',5);                
data_m4.Properties.VariableNames = Headers;

alpha_m4 = data_m4.alpha;
Cm_m4 = data_m4.Cm;
CL_m4 = data_m4.CL; 
CD_m4 = data_m4.CD; 

for i=1:length(alpha_m4)
  alpha_0_m4 = interp1(CL_m4,alpha_m4,0);
end

alpha_zero_torque_m4 = interp1(Cm_m4,alpha_m4-alpha_0_m4,0);

Cl_zero_moment_m4 = interp1(alpha_m4-alpha_0_m4,CL_m4,alpha_zero_torque_m4);

V_m4 = sqrt(2*W/(rho*S*Cl_zero_moment_m4));

%%%%%%% NOTE: The trim velocity is found at approx 13.4 m/s. However, the
%%%%%%% minimum velocity allowed by XFLR5 is 17 m/s in order to avoid 
%%%%%%% convergence errors. Hence, it is assumed 13.4 m/s as the trim velocity
%%%%%%% in this case of study, as convergence has not been reached. 

R_m4 = deltaH.*CL_m4(CL_m4>=0)./CD_m4(CL_m4>=0);

rate_descent_m4 = - sqrt(2*W/(rho*S)).*CD_m4(CL_m4>=0)./(CL_m4(CL_m4>=0).^(3/2)); 

Endurance_m4 = -deltaH./rate_descent_m4;

% Trim condition

CD_zero_moment_m4 = interp1(alpha_m4-alpha_0_m4,CD_m4,alpha_zero_torque_m4);

R_Trim_m4 = deltaH*Cl_zero_moment_m4/CD_zero_moment_m4;

E_Trim_m4 = deltaH/(sqrt(2*W/(rho*S)).*CD_zero_moment_m4/(Cl_zero_moment_m4^(3/2)));

%% Twist +1 

data_1 =readtable('TW1_T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
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

%%%%%%% NOTE: The trim velocity is found at approx 13.4 m/s. However, the
%%%%%%% maximum velocity allowed by XFLR5 is 40 m/s in order to avoid 
%%%%%%% convergence errors. Hence, it is assumed 13.4 m/s as the trim velocity
%%%%%%% in this case of study, as convergence has not been reached. 

R_1 = deltaH.*CL_1(CL_1>=0)./CD_1(CL_1>=0);

rate_descent_1 = - sqrt(2*W/(rho*S)).*CD_1(CL_1>=0)./(CL_1(CL_1>=0).^(3/2)); 

Endurance_1 = -deltaH./rate_descent_1;

% Trim condition

CD_zero_moment_1 = interp1(alpha_1-alpha_0_1,CD_1,alpha_zero_torque_1);

R_Trim_1 = deltaH*Cl_zero_moment_1/CD_zero_moment_1;

E_Trim_1 = deltaH/(sqrt(2*W(1)/(rho*S(1))).*CD_zero_moment_1/(Cl_zero_moment_1^(3/2)));

%% Plots 
figure(1)
hold on 
plot(alpha_Nom-alpha_0_Nom,Cm_Nom)
plot(alpha_m6-alpha_0_m6,Cm_m6)
plot(alpha_m4-alpha_0_m4,Cm_m4)
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
plot(alpha_m6-alpha_0_m6,CL_m6)
plot(alpha_m4-alpha_0_m4,CL_m4)
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
plot(Cm_m6,CL_m6)
plot(Cm_m4,CL_m4)
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
plot(alpha_m6(CL_m6>=0)-alpha_0_m6,Endurance_m6,'r-','LineWidth',1)
plot(alpha_m4(CL_m4>=0)-alpha_0_m4,Endurance_m4,'m-','LineWidth',1)
plot(alpha_1(CL_1>=0)-alpha_0_1,Endurance_1,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,E_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_m6,E_Trim_m6,'r.','MarkerSize',20)
plot(alpha_zero_torque_m4,E_Trim_m4,'m.','MarkerSize',20)
plot(alpha_zero_torque_1,E_Trim_1,'k.','MarkerSize',20)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Endurance [s]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different twist angles at the tip','Interpreter','latex','FontSize',14)

figure(5)
hold on 
plot(alpha_Nom(CL_Nom>=0)-alpha_0_Nom,R_Nom,'b--','LineWidth',1)
plot(alpha_m6(CL_m6>=0)-alpha_0_m6,R_m6,'r-','LineWidth',1)
plot(alpha_m4(CL_m4>=0)-alpha_0_m4,R_m4,'m-','LineWidth',1)
plot(alpha_1(CL_1>=0)-alpha_0_1,R_1,'k-','LineWidth',1)

plot(alpha_zero_torque_Nom,R_Trim_Nom,'b.','MarkerSize',20)
plot(alpha_zero_torque_m6,R_Trim_m6,'r.','MarkerSize',20)
plot(alpha_zero_torque_m4,R_Trim_m4,'m.','MarkerSize',20)
plot(alpha_zero_torque_1,R_Trim_1,'k.','MarkerSize',20)
yline(0,'--')
grid minor
axis square
xlim([0,7.5])
xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex','FontSize',14)
ylabel('$Range [m]$','Interpreter','latex','FontSize',14)
legend(leyenda,'Interpreter','latex')
title('Aircraft performance for different twist angles at the tip','Interpreter','latex','FontSize',14)

figure(6)
hold on
grid minor
axis square
plot([-6,-4,-1,1], [V_m6,V_m4,V_Nom,V_1],'r.-','MarkerSize',10,'LineWidth',2)
xlabel('Twist at the tip [$^{\circ}$]','Interpreter','latex','FontSize',18)
ylabel('Trim velocity [m/s]','Interpreter','latex','FontSize',18)
title('Trim condition for different twist angles at the tip','Interpreter','latex','FontSize',18)
