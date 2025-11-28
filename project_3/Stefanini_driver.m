clear all;
close all;

% Driver: Simplified model for Stefanini et al
% one isoform of VEGF; one receptor; one antibody
% three tissues (blood (central), tumor, rest of body)
% 2025-11-18

% index array (for readability of equations)
n.VEGF_b = 1;
n.VEGFR2_b = 2;
n.VEGFVEGFR2_b = 3;
n.Ab_b = 4;
n.VEGFAb_b = 5;
n.VEGF_t = 6;
n.VEGFR2_t = 7;
n.VEGFVEGFR2_t = 8;
n.Ab_t = 9;
n.VEGFAb_t = 10;
n.VEGF_r = 11;
n.VEGFR2_r = 12;
n.VEGFVEGFR2_r = 13;
n.Ab_r = 14;
n.VEGFAb_r = 15;

% parameters - compartment volumes (specifically, interstitial space)
p.Vol_b = 5*.6; % 5 litres, 60% available
p.Vol_t = 1*.611; % 1 litre, 61% available
p.Vol_r = 40*.0816; % 40 litres, 8% available

% parameters - clearance and internalization
p.kcl_V  = 0.0648; % min-1
p.kcl_VA = 2.2e-5;
p.kcl_A  = 2.2e-5;

p.kintR2  = 0.0168; %min-1
p.kintVR2 = 0.0168;

% parameters - production
p.VEGFprod_b = 0;
p.VEGFprod_t = 3; 
p.VEGFprod_r = 10; 

p.VEGFR2prod_b = 0; 
p.VEGFR2prod_t = 2.85e2*p.kintR2; % p.VEGFR2prod_t = goal conc * p.kintR2
p.VEGFR2prod_r = 2.2e3*p.kintR2; % p.VEGFR2prod_n = goal conc * p.kintR2

% parameters - binding
p.konVR  = 6e-4; % pM-1 min-1 
p.koffVR = 0.06; % min-1

p.konVA  = 5.52e-6; % pM-1 min-1
p.koffVA = 0.012;    % min-1

% compartmental transport
p.k_bt = 4.12e-4; % min-1
p.k_br = 3.19e-3; % min-1
p.k_tb = 4.12e-4; % min-1
p.k_rb = 3.19e-3+6.15e-4; % min-1; additional effect of lymphatic transport
p.AbEx = 1; % set to 0 for no antibody transport; set to 1 for regular transport

% antibody dose
dose = 10*70/150000000*1e12; % 10mg/kg * person size = 70 kg => 700 mg; MW = 150 kDa = 150,000,000 mg/mol

% simulations - 1 - run to steady state, no antibody (observe levels of
% VEGF in each tissue)

% initial conditions 
y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@VEGFAbeqns,[0:1:sstime],y0,options,p,n);
y0ss = Y1(end,:)';

p.AbEx = 0; % set to 0 for no antibody transport; set to 1 for regular transport

y0ss(n.Ab_b)=y0ss(n.Ab_b)+dose/p.Vol_b; % add antibody dose
endtime = 60*24*21; 

[T2,Y2] = ode45(@VEGFAbeqns,[sstime:1:(sstime+endtime)],y0ss,options,p,n);
Tout = [T1;T2]-sstime; % normalize to antibody addition at time zero
Yout = [Y1;Y2];

%% VISUALIZATION

% figure('visible',on); % only show figures if flag is 'on'
f3 = figure;

ax3a=subplot(3,1,1);
plot(ax3a,Tout,Yout(:,n.VEGF_b),'k','linewidth',3)
title(ax3a,'Concentration of free VEGF in Blood')
ylabel(ax3a,'[VEGF] (pM)')
xlabel(ax3a,'time (mins)')

ax3b=subplot(3,1,2);
plot(ax3b,Tout,Yout(:,n.VEGF_t),'k','linewidth',3)
title(ax3b,'Concentration of free VEGF in Tumor')
ylabel(ax3b,'[VEGF] (pM)')
xlabel(ax3b,'time (mins)')

ax3c=subplot(3,1,3);
plot(ax3c,Tout,Yout(:,n.VEGF_r),'k','linewidth',3)
title(ax3c,'Concentration of free VEGF in Rest of Body')
ylabel(ax3c,'[VEGF] (pM)')
xlabel(ax3c,'time (mins)')

f4 = figure;

ax4a=subplot(3,1,1);
plot(ax4a,Tout,Yout(:,n.VEGFR2_b),'k','linewidth',3)
title(ax4a,'Concentration of VEGFR2 in Blood')
ylabel(ax4a,'[VEGFR2] (pM)')
xlabel(ax4a,'time (mins)')

ax4b=subplot(3,1,2);
plot(ax4b,Tout,Yout(:,n.VEGFR2_t),'k','linewidth',3)
title(ax4b,'Concentration of VEGFR2 in Tumor')
ylabel(ax4b,'[VEGFR2] (pM)')
xlabel(ax4b,'time (mins)')

ax4c=subplot(3,1,3);
plot(ax4c,Tout,Yout(:,n.VEGFR2_r),'k','linewidth',3)
title(ax4c,'Concentration of VEGFR2 in Rest of Body')
ylabel(ax4c,'[VEGFR2] (pM)')
xlabel(ax4c,'time (mins)')

f5 = figure;

ax5a=subplot(3,1,1);
plot(ax5a,Tout,Yout(:,n.Ab_b),'k','linewidth',3)
title(ax5a,'Concentration of antibody in Blood')
ylabel(ax5a,'[Ab] (pM)')
xlabel(ax5a,'time (mins)')

ax5b=subplot(3,1,2);
plot(ax5b,Tout,Yout(:,n.Ab_t),'k','linewidth',3)
title(ax5b,'Concentration of antibody in Tumor')
ylabel(ax5b,'[Ab] (pM)')
xlabel(ax5b,'time (mins)')

ax5c=subplot(3,1,3);
plot(ax5c,Tout,Yout(:,n.Ab_r),'k','linewidth',3)
title(ax5c,'Concentration of antibody in Rest of Body')
ylabel(ax5c,'[Ab] (pM)')
xlabel(ax5c,'time (mins)')

f6 = figure;

ax6a=subplot(3,1,1);
plot(ax6a,Tout,Yout(:,n.VEGFAb_b),'k','linewidth',3)
title(ax6a,'Concentration of VEGF-antibody complex in Blood')
ylabel(ax6a,'[Ab-V] (pM)')
xlabel(ax6a,'time (mins)')

ax6b=subplot(3,1,2);
plot(ax6b,Tout,Yout(:,n.VEGFAb_t),'k','linewidth',3)
title(ax6b,'Concentration of VEGF-antibody complex in Tumor')
ylabel(ax6b,'[Ab-V] (pM)')
xlabel(ax6b,'time (mins)')

ax6c=subplot(3,1,3);
plot(ax6c,Tout,Yout(:,n.VEGFAb_r),'k','linewidth',3)
title(ax6c,'Concentration of VEGF-antibody complex in Rest of Body')
ylabel(ax6c,'[Ab-V] (pM)')
xlabel(ax6c,'time (mins)')
