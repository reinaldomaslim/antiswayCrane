%Closed loop crane anti-sway control


%Assumptions:

%Negligible string elongation
%Zero friction
%Small deviation of theta
%Neglect effect of initial jerk 


clc
clear all
close all

global M m l g D Vmax w K T a t1 t2 t_total T1 T2 t e P A1 A2; 

% Default l=40, M=26 tons, m=67.3 tons, Vmax=3m/s


% Crane constants
M = 26*10^3;  % Mass of trolley (kg)
m = 67.3*10^3;  % Mass of load (kg)
l = 40;  % Length of hoist rope (m)
D = 120;  % Desired travel distance (m)

Vmax= 2.8; %Maximum trolley speed (m/s)
g = 9.81;   % Gravitational acceleration(m/s^2)
w = sqrt(g/l); %natural frequency (rad/s)
T= 2*pi/w; %period (s)
K=2/w;
a=0.5; %Constant acceleration (or average)

%ZV filter coeffs

e=0.05; %damping coef
P=exp(-e*pi/sqrt(1-e^2));
A1=1/(1+P);
A2=P/(1+P);
t=pi/(w*sqrt(1-e^2));

t1= round(100*Vmax/a)/100;
t2= round(100*(D-a*t1^2)/Vmax)/100;

T1= t1;
T2= t1+t2;
t_total=T2+t1;


