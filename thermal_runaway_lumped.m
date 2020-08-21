close all
clc
p=2172.99;
c=1389.7;
K=27.5;
h1=257;
h2=1714;
h3=81;
h4=155;
w1=6.1*10^5;
w2=6.1*10^5;
w3=1.221*10^6;
w4=4.069*10^5;
a1=1.667*10^15;
a2=2.5*10^13;
a3=6.667*10^13;
a4=5.14*10^25;
% a1=(2.57*10^5)*(1.667*10^15)*(1.39*10^3);
% a2=(1.714*10^6)*(2.5*10^13)*(1.39*10^3);
% a3=(3.14*10^5)*(6.667*10^13)*(1.3*10^3);
% a4=(1.55*10^5)*(5.14*10^25)*(500);
e1=1.3508*10^5;
e2=1.3508*10^5;
e3= 1.396*10^5;
e4=2.74*10^5;
R=8.314;
T0=293;
dt=1;
t=5000;
X=zeros(1,t);
T=zeros(1,t);
S=zeros(1,t);
F=zeros(1,t);
DT=T;
S(1)=T0;
T(1)=T0;
I=80;
r=0.038;
sigma=5.6*10^(-8);
c1=0.15;
c2=0.75;
c3=0.04;
c4=1;

for i=1:1:t
    if T(i)<378
        q=0;
    elseif T(i)>=378&&T(i)<393
        k1=(a1*c1*exp(-e1/(R*T(i))));
        q=w1*k1*h1;
        c1=c1+(-k1)*dt;
        if c1<0
            c1=0;
        end
    elseif T(i)>=393&&T(i)<443
        k1=(a1*c1*exp(-e1/(R*T(i))));
        q1=w1*k1*h1;
        k2=(a2*c2*exp(-e2/(R*T(i))));
        q2=w2*k2*h2;
        q=q1+q2;
        c1=c1-k1*dt;
        c2=c2-k2*dt;
        if c1<0
            c1=0;
        end
        if c2<0
            c2=0;
        end
    elseif T(i)>=443&&T(i)<473
        k1=(a1*c1*exp(-e1/(R*T(i))));
        q1=w1*k1*h1;
        k2=(a2*c2*exp(-e2/(R*T(i))));
        q2=w2*k2*h2;
        k3=(a3*c3*(1-c3)*exp(-e3/(R*T(i))));
        q3=w3*k3*h3;
        q=q1+q2+q3;
        c1=c1-k1*dt;
        c2=c2-k2*dt;
        c3=c3-k3*dt;
        if c1<0
            c1=0;
        end
        if c2<0
            c2=0;
        end
        if c3<0
            c3=0;
        end
    else
        k1=(a1*c1*exp(-e1/(R*T(i))));
        q1=w1*k1*h1;
        k2=(a2*c2*exp(-e2/(R*T(i))));
        q2=w2*k2*h2;
        k3=(a3*c3*(1-c3)*exp(-e3/(R*T(i))));
        q3=w3*k3*h3;
        k4=(a4*c4*exp(-e4/(R*T(i))));
        q4=w4*k4*h4;
        q=q1+q2+q3+q4-(dt/(p*c))*(((I^2)*r)/(pi*81*4225*10^(-9)));
        c1=c1-k1*dt;
        c2=c2-k2*dt;
        c3=c3-k3*dt;
        c4=c4-k4*dt;
        if c1<0
            c1=0;
        end
        if c2<0
            c2=0;
        end
        if c3<0
            c3=0;
        end
        if c4<0
            c4=0;
        end

    end
    T(i+1)=T(i)+(dt/(p*c))*(((I^2)*r)/(pi*81*4225*10^(-9))+q-((sigma)*(T(i)^4-T0^4)));
    DT(i)=(T(i+1)-T(i))/dt;
    
     if S(i)<378
        q1=0;
    elseif S(i)>=378&&S(i)<393
        q1=(a1*exp(-e1/(R*S(i))));
    elseif S(i)>=393&&S(i)<443
        q1=(a2*exp(-e2/(R*S(i))));
    elseif S(i)>=443&&S(i)<473
        q1=(a3*exp(-e3/(R*S(i))));
    else
        q1=(a4*exp(-e4/(R*S(i)))); 
    end
    S(i+1)=S(i)+(dt*K/(p*c*0.01))*(T(i)-S(i))+((sigma/(c*p*0.01))*(T(i)^4-T0^4));
    X(i)=i;
    
end
for i=2:1:t
F(i)=(T(i+1)-2*T(i)+T(i-1))/dt^2;
end
F(1)=0;
X(t+1)=t+1;
DT(t+1)=(T(t+1)-T(t))/dt;
figure(1)
plot(X,T)
hold on
plot(X,S)
xlabel('time')
ylabel('temperature of cells')
% ylim([280 1000])
title('T_t graph')
legend('first cell','second cell')

figure(2)
 plot(T,DT)
 xlabel('temperature')
ylabel('DT/dt')
title('DT_T graph')
 

 
 
  


    
    

    

