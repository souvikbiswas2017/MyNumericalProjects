clear;clc;
%one-dimensional euler equations with a shock
%grid and time step
nxp=101;
nstep=5000;
length=100.0;
h=length/(nxp-1);
alpha=10.0;
dt=0.00005;
%output varibles
x=zeros(nxp,1);
rho=zeros(nxp,1);
u=zeros(nxp,1);
p=zeros(nxp,1);
e=zeros(nxp,1);
M=zeros(nxp,1);
E=zeros(nxp,1);
H=zeros(nxp,1);
P=zeros(nxp,1);
Ma=zeros(nxp,1);

utmp=zeros(nxp-1,1);
ptmp=zeros(nxp-1,1);
rhotmp=zeros(nxp-1,1);
Mtmp=zeros(nxp-1,1);
Etmp=zeros(nxp-1,1);

time=0.0;
%system parameters
gamma=1.4;
%initial left and right state
rhol=10;
rhor=1;
pl=1e6;
pr=1e5;
ul=100.0;
ur=0.0;
%problem set up and initial conditions
for i=2:nxp
    x(i)=(i-1)*h;
end;

rho(1)= rhol;
p(1) = pl;
e(1)=p(1)/rho(1)/(gamma-1);
u(1)=ul;
for i=2:nxp
    rho(i)=rhor;
    p(i)=pr;
    e(i)=p(i)/rho(i)/(gamma-1);
    u(i)=ur;
end;
% while(time<finaltime)
for cnt=1:nstep
    for i=1:nxp
        M(i)=rho(i)*u(i);
        E(i)=rho(i)*e(i)+0.5*rho(i)*u(i)^2;
        H(i)=p(i)+rho(i)*u(i)^2;
        P(i)=u(i)*(rho(i)*(e(i)+0.5*u(i)^2)+p(i));
        Ma(i)=u(i)/sqrt(gamma*p(i)/rho(i));
    end;
	for i=1:nxp-1
        rhotmp(i)=0.5*(rho(i)+rho(i+1))-0.5*(dt/h)*(M(i+1)-M(i));
        Mtmp(i)=0.5*(M(i)+M(i+1))-0.5*(dt/h)*(H(i+1)-H(i));
        Etmp(i)=0.5*(E(i)+E(i+1))-0.5*(dt/h)*(P(i+1)-P(i));

        utmp(i)=Mtmp(i)/rhotmp(i);
        ptmp(i)=(gamma-1)*rhotmp(i)*(Etmp(i)/rhotmp(i)-0.5*utmp(i)^2);
    end;
    for i=2:nxp-1
        rho(i)=rho(i)-(dt/h)*(rhotmp(i)*utmp(i)-rhotmp(i-1)*utmp(i-1));
        M(i)=M(i)-(dt/h)*(Mtmp(i)*utmp(i)+ptmp(i)-Mtmp(i-1)*utmp(i-1)-ptmp(i-1));
        M(i)=M(i)+alpha*(dt/h)*(rho(i)*abs(u(i+1)-u(i))*(u(i+1)-u(i))-...
                                rho(i-1)*abs(u(i)-u(i-1))*(u(i)-u(i-1)));
        E(i)=E(i)-(dt/h)*(Etmp(i)*utmp(i)+ptmp(i)*utmp(i)-Etmp(i-1)*utmp(i-1)-ptmp(i-1)*utmp(i-1));
        E(i)=E(i)+alpha*(dt/h)*(rho(i)*u(i)*abs(u(i+1)-u(i))*(u(i+1)-u(i))-...
                                rho(i-1)*u(i-1)*abs(u(i)-u(i-1))*(u(i)-u(i-1)));
        u(i)=M(i)/rho(i);
        e(i)=(E(i)/rho(i))-0.5*u(i)^2;
        p(i)=(gamma-1)*rho(i)*e(i);
        Ma(i)=u(i)/sqrt(gamma*p(i)/rho(i));
    end;
    
    for i=2:nxp-1
        rho(i) = 0.5*(rho(i+1)+rho(i-1));
    end;
    time=time+dt;

    subplot(221);
	plot(x,rho);
	title(sprintf('Density Plot at t = %5.3f for %dX1 resolution',time,nxp-1));
	xlabel('X');ylabel('Density');
    %axis([0.0 10.0 0.0 1.25]);
    
    subplot(222);
	plot(x,u);
	title(sprintf('Velocity Plot at t = %5.3f for %dX1 resolution',time,nxp-1));
	xlabel('X');ylabel('Velocity');
    %axis([0.0 10.0 0.0 400.0]);

    subplot(223);
	plot(x,p);
	title(sprintf('Pressure Plot at t = %5.3f for %dX1 resolution',time,nxp-1));
	xlabel('X');ylabel('Pressure');
    %axis([0.0 10.0 0.0 1.25e5]);

    subplot(224);
	plot(x,Ma);
	title(sprintf('Mach No. Plot at t = %5.3f for %dX1 resolution',time,nxp-1));
	xlabel('X');ylabel('Ma');
    %axis([0.0 10.0 0.0 1.25]);

    pause(0.001);
end;
    