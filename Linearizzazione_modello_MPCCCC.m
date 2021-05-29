clear all
clc

% Declare symbolic variables
syms n mu vx vy r Fm delta dFm dDelta param_lf param_lr param_m Mz param_Iz tire_forces_rearF_y tire_forces_frontF_y param_Cr0 param_Cr2 k
syms param_Df param_Cf param_Bf alpha_f alpha_r param_Dr param_Cr param_Br dk

alpha_f = sym("alpha_f","real");
alpha_r = sym("alpha_r","real");
vx = sym("vx","real");
vy = sym("vy","real");
param_lf = sym("param_lf","real");
r = sym("r","real");


% Friction force
friction_force = - param_Cr0 - param_Cr2*vx*vx;

% Tire front force
alpha_f = -atan2(vy+r*param_lf,vx) + delta;
tire_forces_frontF_y = param_Df * sin(param_Cf * atan(param_Bf * alpha_f ));

% Tire rear force
tire_forces_rearF_y = param_Dr * sin(param_Cr * atan(param_Br * alpha_r ));
alpha_r = atan2(vy-r*param_lr,vx);

% Continuous time model
f1 = vx * sin(mu) + vy * cos(mu);
f2 = r - k * (vx * cos(mu) - vy * sin (mu)/(1-n*k));
f3 = 1.0/param_m*(Fm * (1+cos(delta)) - friction_force - tire_forces_frontF_y*sin(delta) + param_m*vy*r);
f4 = 1.0/param_m*(Fm * sin(delta) + tire_forces_rearF_y + tire_forces_frontF_y*cos(delta) - param_m*vx*r);
f5 = 1.0/param_Iz*(Fm * sin(delta)*param_lf + tire_forces_frontF_y*param_lf*cos(delta) - tire_forces_rearF_y*param_lr + Mz);
f6 = dFm;
f7 = dDelta;
f8 = dk;

% Matrix A e B

f = f2;

df = ccode(diff(f,n));
df = ccode(diff(f,mu));
df = ccode(diff(f,vx));
df = ccode(diff(f,vy));
df = ccode(diff(f,r));
df = ccode(diff(f,Fm));
df = ccode(diff(f,delta));
df = ccode(diff(f,k))

df = ccode(diff(f,dFm));
df = ccode(diff(f,dDelta));
df = ccode(diff(f,Mz));
df = ccode(diff(f,dk))










