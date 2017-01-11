%Open Loop Indoor crane anti-sway control


%Assumptions:
%Zero disturbances
%Negligible string elongation
%Zero friction
%Small deviation of theta
%Neglect effect of initial jerk 


clc
clear all
close all

global M m l g D t1 t2 Vmax w T thetamax a_max a n t_total V_average T1 T2; 

% Default l=40, M=26 tons, m=67.3 tons, Vmax=3m/s


% Crane constants
M = 26*10^3;  % Mass of trolley (kg)
m = 67.3*10^3;  % Mass of load (kg)
l = 40;  % Length of hoist rope (m)
D = 100;  % Desired travel distance (m)

Vmax= 3; %Maximum trolley speed (m/s)
g = 9.81;   % Gravitational acceleration(m/s^2)
w = sqrt(g/l); %natural frequency (rad/s)
T= 2*pi/w; %period (s)
thetamax= 5*(pi/180); %max sway angle(rad)

%Velocity profile parameters(s)
a_max=thetamax*g/2; %Estimated max a for thetamax=5deg (m/s^2)
t1 = T/4;% pulse width/non-pulse width (s)

%Velocity profile parameters for t2>0(s)

n = ceil(Vmax/(2*a_max*t1)); %number of +ve/-ve pulse pairs, take ceiling to reach Vmax
a = Vmax/(2*n*t1); %Corresponding acceleration for n
t2 = (D-8*(n^2)*(t1^2)*a)/(2*n*a*t1);  % width when trolley travelling at Vmax(s)

if t2<0
    a=D/(8*(t1^2)*(n^2));  %Corresponding a for desired distance D 
    t2=0;
end 




t_total = 2*n*T + t2; %Total time for travel
V_average = D/t_total; %Average velocity

T1= n*T; %Alternating pulses time
T2= n*T+t2; % before deccelerate time


