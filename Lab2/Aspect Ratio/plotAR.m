clear; clc; close all; 
%% Global parameters 

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

leyenda = {'Nominal Case','AR = 6.89','AR = 13.79','AR = 27.59'};


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

R_Nom = deltaH.*CL_Nom(CL_Nom>=0)./CD_Nom(CL_Nom>=0);

rate_descent_Nom = - sqrt(2*W(1)/(rho*S(1))).*CD_Nom(CL_Nom>=0)./(CL_Nom(CL_Nom>=0).^(3/2)); 

Endurance_Nom = -deltaH./rate_descent_Nom;

%% AR = 6.89 
data_6 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_6.Properties.VariableNames = Headers;

alpha_6 = data_6.alpha;
Cm_6 = data_6.Cm;
CL_6 = data_6.CL; 
CD_6 = data_6.CD; 

for i=1:length(alpha_6)
  alpha_0_6 = interp1(CL_6,alpha_6,0);
end

R_6 = deltaH.*CL_6(CL_6>=0)./CD_6(CL_6>=0);

rate_descent_6 = - sqrt(2*W(2)/(rho*S(2))).*CD_6(CL_6>=0)./(CL_6(CL_6>=0).^(3/2)); 

Endurance_6 = -deltaH./rate_descent_6;

%% AR = 13.793
data_13 =readtable('AR_13_793_T1-20_6 m_s-VLM2.txt','HeaderLines',5);                
data_13.Properties.VariableNames = Headers;

alpha_13 = data_13.alpha;
Cm_13 = data_13.Cm;
CL_13 = data_13.CL; 

for i=1:length(alpha_13)
  alpha_0_13 = interp1(CL_13,alpha_13,0);
end

%% AR = 13.89
data_13 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_13.Properties.VariableNames = Headers;

alpha_13 = data_13.alpha;
Cm_13 = data_13.Cm;
CL_13 = data_13.CL; 
CD_13 = data_13.CD; 

for i=1:length(alpha_13)
  alpha_0_13 = interp1(CL_13,alpha_13,0);
end

R_13 = deltaH.*CL_13(CL_13>=0)./CD_13(CL_13>=0);

rate_descent_13 = - sqrt(2*W(3)/(rho*S(3))).*CD_13(CL_13>=0)./(CL_13(CL_13>=0).^(3/2)); 

Endurance_13 = -deltaH./rate_descent_13;

%% AR = 27.59
data_27 =readtable('T1-20_0 m_s-VLM2.txt','HeaderLines',5);                
data_27.Properties.VariableNames = Headers;

alpha_27 = data_27.alpha;
Cm_27 = data_27.Cm;
CL_27 = data_27.CL; 
CD_27 = data_27.CD; 

for i=1:length(alpha_27)
  alpha_0_27 = interp1(CL_27,alpha_27,0);
end

R_27 = deltaH.*CL_27(CL_27>=0)./CD_27(CL_27>=0);

rate_descent_27 = - sqrt(2*W(1)/(rho*S(1))).*CD_27(CL_27>=0)./(CL_27(CL_27>=0).^(3/2)); 

Endurance_27 = -deltaH./rate_descent_27;


%% Plots 
% figure(1)
% hold on 
% plot(alpha-alpha_0,Cm)
% yline(0,'--')
% grid minor
% axis square
% xlim([-4,7.5])
% xlabel('$\alpha$ [rad]','Interpreter','latex')
% ylabel('$C_{m}$','Interpreter','latex')
% legend(leyenda,'Interpreter','latex')
% title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')
% 
% figure(2)
% hold on 
% plot(alpha-alpha_0,CL)
% yline(0,'--')
% grid minor
% axis square
% xlim([-4,7.5])
% xlabel('$\alpha$ [rad]','Interpreter','latex')
% ylabel('$C_{L}$','Interpreter','latex')
% legend(leyenda,'Interpreter','latex')
% title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')
% 
% figure(3)
% hold on 
% plot(CL,Cm)
% yline(0,'--')
% xlim([-0.5,0.6])
% grid minor
% axis square
% xlabel('$C_{L}$' ,'Interpreter','latex')
% ylabel('$C_{m}$','Interpreter','latex')
% legend(leyenda,'Interpreter','latex')
% title('Analysis of the aircraft performance for different $x_{CG}$','Interpreter','latex')
