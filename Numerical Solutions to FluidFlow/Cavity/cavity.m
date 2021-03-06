
%Two Dimensional Driven Cavity Problem by MAC Method
clear;clc;
%Grid variables
nx=64;ny=64;
xlength=1;ylength=1;
h=xlength/nx;
dt=0.0003;nstep=50;
%Problem parameters
uwall=1;rho=1;meu=0.1;
%General variables
beta=1.5;epsP=0.0;eps=1.0;
%Initializing the Arrays
P=zeros(nx+2,ny+2);
u=zeros(nx+1,ny+2);v=zeros(nx+2,ny+1);
ut=zeros(nx+1,ny+2);vt=zeros(nx+2,ny+1);
uu=zeros(nx+1,nx+1);vv=zeros(nx+1,nx+1);pp=zeros(nx,ny);w=zeros(nx+1,nx+1);
p0=zeros(nx+2,nx+2);a=zeros(nx+2,ny+2);sf=zeros(nx+1,ny+1);
for i=1:nx+2
    for j=1:ny+2
        if (i==1 | j==1 | i==nx+2 | j==ny+2)
            a(i,j)=0.0;
        elseif (i==2 | i==nx+1) & (j==2 | j==ny+1)
            a(i,j)=1/2.0;
        elseif (i==2 | i==nx+1 | j==2 | j==ny+1)
            a(i,j)=1/3.0;
        else
            a(i,j)=1/4.0;
        end;
    end;
end;
%Time Loop
for k=1:nstep
    %Velocity Boundary Conditions
    u(1:nx+1,ny+2)=2*uwall-u(1:nx+1,ny+1);
    u(1:nx+1,1)=-u(1:nx+1,2);
    v(nx+2,1:ny+1)=-v(nx+1,1:ny+1);
    v(1,1:ny+1)=-v(2,1:ny+1);
    %Finding Projected Velocity in X
    for i=2:nx
        for j=2:ny+1
            ut(i,j)=u(i,j)+dt*(-(0.25/h)*((u(i+1,j)+u(i,j))^2-(u(i,j)+...
                    u(i-1,j))^2+(u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+...
                    (meu/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)));
        end;
    end;
    %Finding Projected Velocity in Y
    for i=2:nx+1
        for j=2:ny
            vt(i,j)=v(i,j)+dt*(-(0.25/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+...
                    (v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)+...
                    (meu/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)));
        end;
    end;
    %Solving the Pressure Equation
    while eps>1e-3
        p0(1:nx+2,1:ny+2)=P(1:nx+2,1:ny+2);
        for i=2:nx+1
            for j=2:ny+1
                P(i,j)=beta*a(i,j)*(P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1)-...
                       (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+(1-beta)*P(i,j);
            end;
        end;
        eps=0;
        for i=2:nx+1
            for j=2:ny+1
                eps = eps+(P(i,j)-p0(i,j))^2;
            end;
        end;
        eps=sqrt(eps/(nx)*(ny));
    end;
    eps=1;
    %Updating the Velocity
    for i=2:nx
        for j=2:ny+1
            u(i,j) = ut(i,j)-(dt/h)*(P(i+1,j)-P(i,j));
        end;
    end;
    for i=2:nx+1
        for j=2:ny
            v(i,j) = vt(i,j)-(dt/h)*(P(i,j+1)-P(i,j));
        end;
    end;
    
    %Calculating vorticity at the pressure points
    for i=1:nx+1
        for j=1:ny+1
            uu(i,j)=0.5*(u(i,j+1)+u(i,j));
            vv(i,j)=0.5*(v(i+1,j)+v(i,j));
            w(i,j)=(1/2*h)*(u(i,j+1)-u(i,j)-v(i+1,j)+v(i,j));
        end;
    end;
    pp(1:nx,1:ny) = P(2:nx+1,2:ny+1);
    for i=2:nx
        for j=2:ny
            sf(i,j)=-0.5*(u(i,j)+u(i-1,j))*h+0.5*(v(i,j)+v(i,j-1))*h+sf(i-1,j-1);
        end;
    end;
    cavityplot;
    k
end;
%calculating 2nd order norm of velocity vector at 9 X 9 at time =0.15
%grid points using constant dt to prove spatial convergence
unorm=0.0;vnorm=0.0;
incx=nx/8;incy=ny/8;
for i=1:incx:nx+1
    for j=1:incy:ny+1
        unorm=unorm+uu(i,j)^2;
        vnorm=vnorm+vv(i,j)^2;
    end;
end;
unorm=incx*incy*sqrt(unorm)/nx/ny;
vnorm=incx*incy*sqrt(vnorm)/nx/ny;
