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

global M m l g D t1 t2 Vmax w T thetamax a_max a a1 n tmid t_total t_total1 t_total2 V_average V2 T1 T2 k p t3 a2; 

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

%Max velocity solution (params: k, tmid, a1, t_total1)

k = ceil(Vmax/(2*a_max*t1)); %number of +ve/-ve pulse pairs, take ceiling to reach Vmax
a1 = Vmax/(2*k*t1); %Corresponding acceleration for n
tmid = (D-(4*k+1)*t1*Vmax)/Vmax;  % width when trolley travelling at Vmax(s)

if tmid<0
    a1=D/(2*(4*k+1)*k*(t1^2));  %Corresponding a for desired distance D 
    tmid=0;
end 

t_total1 = 2*k*T + tmid; %Total time for travel

%Max acceleration solution (params: p, t3, a2, t_total2)

p=floor(Vmax/(2*a_max*t1));
V2=2*p*t1*a_max;
a2=a_max;
t3 = (D-(4*p+1)*t1*V2)/V2;


if t3<0
    a2=D/(2*(4*p+1)*p*(t1^2));  %Corresponding a for desired distance D 
    t3=0;
end 

t_total2 = 2*p*T + t3; %Total time for travel

%Comparing both time optimization alternatives
if t_total1<t_total2
    n=k;
    a=a1;
    t_total=t_total1;
    t2=tmid;
else
    n=p;
    a=a2;
    t_total=t_total2;
    t2=t3;
end


V_average = D/t_total; %Average velocity

T1= n*T; %Alternating pulses time
T2= n*T+t2; % before deccelerate time


