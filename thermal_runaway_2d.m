close all
clc
%% Properties
rho=2172.99;
cp=1389.7;
kr=27.5;
ky=0.898;
radius=0.009;
height=65*10^(-3);

%% Reaction constansts
h1=257;
h2=1714;
h3=625;
h4=155;

w1=6.1*10^5;
w2=6.1*10^5;
w3=1.221*10^6;
w4=4.069*10^5;

a1=1.667*10^15;
a2=2.5*10^13;
a3=6.667*10^13;
a4=5.14*10^25;

e1=1.3508*10^5;
e2=1.3508*10^5;
e3= 1.396*10^5;
e4=2.74*10^5;

c1=0.15;
c2=0.75;
c3=0.04;
c4=1;
%% Heat transfer constants
h=20;
sig=5.67*10^-8;
e=0.5;

%% Constants
To=298;
R_univ=8.314;
r=1;
i=80;
area=pi*radius^2*height;
i_eng=(r*i^2)/area;

%% Intialization
dt=0.001;
t=2;
nt=t/dt;
dx=0.001;
lx=radius;
dy=0.001;
ly=height;
nx=lx/dx+1;
ny=ly/dy+1;

time=zeros(1,nt);
radial=zeros(1,nx);

T=zeros(nt,nx,ny);
%% Initial Conditions
for i=1:1:nx
    for j=1:1:ny
        T(1,i,j)=To;   
    end
end
%% Computations

for i=1:nt
            % energy generation
            if T(i,1,1)<378
                q=0;
            elseif T(i,1,1)>=378&&T(i,1,1)<393
                k1=(a1*c1*exp(-e1/(R*T(i,1,1))));
                q=w1*k1*h1;
                c1=c1+(-k1)*dt;
                if c1<0
                    c1=0;
                end
            elseif T(i,1,1)>=393&&T(i,1,1)<443
                 k1=(a1*c1*exp(-e1/(R*T(i,1,1))));
                q1=w1*k1*h1;
                k2=(a2*c2*exp(-e2/(R*T(i,1,1))));
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
            elseif T(i,1,1)>=443&&T(i,1,1)<473
                k1=(a1*c1*exp(-e1/(R*T(i,1,1))));
                q1=w1*k1*h1;
                k2=(a2*c2*exp(-e2/(R*T(i,1,1))));
                q2=w2*k2*h2;
                k3=(a3*c3*(1-c3)*exp(-e3/(R*T(i,1,1))));
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
                k1=(a1*c1*exp(-e1/(R*T(i,1,1))));
                q1=w1*k1*h1;
                k2=(a2*c2*exp(-e2/(R*T(i,1,1))));
                q2=w2*k2*h2;
                k3=(a3*c3*(1-c3)*exp(-e3/(R*T(i,1,1))));
                q3=w3*k3*h3;
                k4=(a4*c4*exp(-e4/(R*T(i,1,1))));
                q4=w4*k4*h4;
                q=q1+q2+q3+q4;
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
     for j=2:nx-1
        for v=2:ny-1
            syms w
            %Interior points
            eqn= (rho*cp)*((w-T(i,j,v))/dt)==(kr/ly)*((T(i,j+1,v)-2*T(i,j,v)+T(i,j-1,v))/(dx)^2)+(ky/lx)*((T(i,j,v+1)-2*T(i,j,v)+T(i,j,v-1))/(dy)^2);
            T(i+1,j,v)=solve(eqn,w);
        end
     end
            % Right Boundary
     for v=1:ny
            if j==nx
                eqn=-k*((w-T(i+1,j-1,v))/dx)==h*(w-To)+sig*e*(w^4-To^4);
                T(i+1,j,v)=solve(eqn,w);
            end
            %Left Boundary
            if j==1
                eqn=-k*((T(i+1,j+1,v)-w)/dx)==q/height+i_eng;
                T(i+1,j,v)=solve(eqn,w);
            end
     end
end
plot(time',T(:,1,ny/2));
