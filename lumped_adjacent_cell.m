close all
clc
m=50*10^(-3);
c=1188.952;
k=0.02;
a1=(2.57*10^5)*(1.667*10^15)*(1.39*10^3);
a2=(1.714*10^6)*(2.5*10^13)*(1.39*10^3);
a3=(3.14*10^5)*(6.667*10^13)*(1.3*10^3);
a4=(1.55*10^5)*(5.14*10^25)*(500);
e1=1.3508*10^5;
e2=1.3508*10^5;
e3= 1.396*10^5;
e4=2.74*10^5;
R=8.314;
T0=378;
dt=0.01;
t=3.5;
I=12;
r=0.038;
vol=pi*81*4225*10^(-9);
a=pi*(0.009)*(0.065);
sigma=5.6*10^(-8);
l=0.002;
T=zeros(t/dt,1);
S=T;
S(1)=311;
T(1)=T0;
X=zeros(t/dt,1);
for i=1:1:t/dt-1
 q=(a1*exp(-e1/(R*T(i))))+ (a2*exp(-e2/(R*T(i)))) +(a3*exp(-e3/(R*T(i)))) + (a4*exp(-e4/(R*T(i)))); 
 T(i+1)=T(i)+(dt/(m*c))*(((I^2)*r)+q*(vol)); 
S(i+1)=S(i)+(dt*a/(m*c))*(((k/l)*(T(i)-S(i)))+sigma*((T(i))^4-(S(i))^4));
%S(i+1)=S(i)+(dt*k*(pi*0.065*0.009)/(m*c*l))*(T(i)-S(i))+((dt*sigma*(pi*0.065*0.009)/(c*m))*(T(i)^4-S(i)^4));
 X(i)=i*dt;
end
X(t/dt)=t;
figure(1)
plot(X,T);
ylabel('Temperature of cell undergoing thermal runaway');
xlabel('Time after start of exothermic reactions');
title('Temperature vs Time graph');
figure(2)
plot(X,S);
ylabel('Temperature of the adjacent cell');
xlabel('Time after start of exothermic reactions');
title('Temperature vs Time graph');



 

 
 
  


    
    

    


