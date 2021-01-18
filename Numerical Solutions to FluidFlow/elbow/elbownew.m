%resolution & parameters governing the run
clear;clc;MaxStep=50; alpha=0.1; dt=0.001;         
%parameters for SOR iteration and grid
MaxIt=1000; Beta=1.5; MaxErr=0.001;                   
h=1.0;t=0.0;Q=1.0;

%-----------------------------------------------------------
%Start Grid Generation
%-----------------------------------------------------------
nB=2;%Number of blocks
%lower boundary
xp(1,1:nB+1)=[0,0,1.5];yp(1,1:nB+1)=[1,0,0];
%upper boundary
xp(2,1:nB+1)=[0.5,0.5,1.5];yp(2,1:nB+1)=[1,0.5,0.5];

Nblock(1:nB)=[8,8];NTot=0;
Mblock(1:nB)=[6,6];MTot=6;

for cnt=1:nB
  N=Nblock(cnt);M=Mblock(cnt);
  x0=xp(1,cnt);x1=xp(1,cnt+1);x3=xp(2,cnt);x2=xp(2,cnt+1);
  y0=yp(1,cnt);y1=yp(1,cnt+1);y3=yp(2,cnt);y2=yp(2,cnt+1);
  Nf=2;if(cnt==1),Nf=1;end
  for i=Nf:N, for j=1:M	
      ii=i+NTot;
      x(ii,j)=((M-j)/(M-1))*(((N-i)/(N-1))*x0+((i-1)/(N-1))*x1)+...
             ((j-1)/(M-1))*(((N-i)/(N-1))*x3+((i-1)/(N-1))*x2);
      y(ii,j)=((M-j)/(M-1))*(((N-i)/(N-1))*y0+((i-1)/(N-1))*y1)+...
             ((j-1)/(M-1))*(((N-i)/(N-1))*y3+((i-1)/(N-1))*y2);
    end, end;
  NTot=NTot+Nblock(cnt)-1;
end
NTot=NTot+1;
%Plot grid
% figure(1);
% axis('equal'), hold on
% for j=1:MTot;plot(x(1:NTot,j),y(1:NTot,j),'linewidth',2);end;
% for i=1:NTot;plot(x(i,1:MTot),y(i,1:MTot),'linewidth',2);end;
% axis('square'),axis([-0.10, 2, 0, 2]); hold on
% axis off
%Smooth grid
for m=1:10
	for i=2:NTot-1
	for j=2:MTot-1
		x(i,j)=0.25*(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1));
		y(i,j)=0.25*(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1));
	end;
	end;
end;
% axis('square'),axis([-0.10, 2, 0, 2]); hold on
% for j=1:MTot;plot(x(1:NTot,j),y(1:NTot,j),'r','linewidth',2);end;
% for i=1:NTot;plot(x(i,1:MTot),y(i,1:MTot),'r','linewidth',2);end;
% axis off
%-----------------------------------------------------------
%End Grid Generation
%-----------------------------------------------------------

nz=NTot;ne=MTot;
sf=zeros(nz,ne);vt=zeros(nz,ne);w=zeros(nz,ne);

%calculate dxde,dxdz,dyde,dydz
for j=1:ne
    for i=2:nz-1
        xz(i,j)=(x(i+1,j)-x(i-1,j))/2/h;
        yz(i,j)=(y(i+1,j)-y(i-1,j))/2/h;
    end;
    xz(1,j)=(x(2,j)-x(1,j))/h;    
    yz(1,j)=(y(2,j)-y(1,j))/h;    
    xz(nz,j)=(x(nz,j)-x(nz-1,j))/h;    
    yz(nz,j)=(y(nz,j)-y(nz-1,j))/h;    
end;
for i=1:nz
    for j=2:ne-1
        xe(i,j)=(x(i,j+1)-x(i,j-1))/2*h;
        ye(i,j)=(y(i,j+1)-y(i,j-1))/2*h;
    end;
    xe(i,1)=(x(i,2)-x(i,1))/h;    
    ye(i,1)=(y(i,2)-y(i,1))/h;    
    xe(i,ne)=(x(i,ne)-x(i,ne-1))/h;    
    ye(i,ne)=(y(i,ne)-y(i,ne-1))/h;
end;

%calculate q1,q2,q3,Jacobian
for i=1:nz
    for j=1:ne
        q1(i,j)=xe(i,j)^2+ye(i,j)^2;
        q2(i,j)=xe(i,j)*xz(i,j)+ye(i,j)*yz(i,j);
        q3(i,j)=xz(i,j)^2+yz(i,j)^2;
        J(i,j)=xz(i,j)*ye(i,j)-xe(i,j)*yz(i,j);
    end;
end;
%calculate the 2nd order derivatives
%xzz,xee,xze,yzz,yez,yee
for i=2:nz-1
    for j=2:ne-1
        xzz(i,j)=(xz(i+1,j)-xz(i-1,j))/2*h;
        yzz(i,j)=(yz(i+1,j)-yz(i-1,j))/2*h;
        yez(i,j)=(ye(i+1,j)-ye(i-1,j))/2*h;
        xez(i,j)=(xe(i+1,j)-xe(i-1,j))/2*h;
        
        xee(i,j)=(xe(i,j+1)-xe(i,j-1))/2*h;
        yee(i,j)=(ye(i,j+1)-ye(i,j-1))/2*h;
        xze(i,j)=(xz(i,j+1)-xz(i,j-1))/2*h;
        yze(i,j)=(yz(i,j+1)-yz(i,j-1))/2*h;
    end;
end;

%calculate laplacian of z and e
for i=2:nz-1
    for j=2:ne-1
        Lz(i,j)=(q1(i,j)*(xe(i,j)*yzz(i,j)-ye(i,j)*xzz(i,j))-...
                2*q2(i,j)*(xe(i,j)*yze(i,j)-ye(i,j)*xze(i,j))+...
                q3(i,j)*(xe(i,j)*yee(i,j)-ye(i,j)*xee(i,j)))/J(i,j)^3;
            
        Le(i,j)=(q1(i,j)*(yz(i,j)*xzz(i,j)-xz(i,j)*yzz(i,j))-...
                2*q2(i,j)*(yz(i,j)*xze(i,j)-xz(i,j)*yze(i,j))+...
                q3(i,j)*(yz(i,j)*xee(i,j)-xz(i,j)*yee(i,j)))/J(i,j)^3;
    end;
end;

for istep=1:MaxStep
    istep
    % define boundary conditons for stream function
    sf(1:nz,1)=0;
    sf(1:nz,ne)=Q;
    for j=2:ne-1
        sf(1,j)=Q*(3*((j-1)/(ne-1))^2-2*((j-1)/(ne-1))^3);
        sf(nz,j)=sf(nz-1,j);
    end;
    
    %SOR solution of stream function
    for iter=1:MaxIt
        sfold=sf;
        for i=2:nz-1
            for j=2:ne-1
                term1 = q2(i,j)*(sf(i+1,j+1)-sf(i-1,j+1)-sf(i+1,j-1)+sf(i-1,j-1))/4/h^2/J(i,j)^2-...
                        q1(i,j)*(sf(i+1,j)+sf(i-1,j))/h^2/J(i,j)^2-...
                        q3(i,j)*(sf(i,j+1)+sf(i,j-1))/h^2/J(i,j)^2-...
                        Lz(i,j)*(sf(i+1,j)-sf(i-1,j))/2/h-...
                        Le(i,j)*(sf(i,j+1)-sf(i,j-1))/2/h-...
                        vt(i,j);
                term2 = -2*(q1(i,j)+q3(i,j))/h^2/J(i,j)^2;                
                sf(i,j) = Beta*(term1/term2)+(1-Beta)*sf(i,j);
            end;
        end;
        Err=0.0;
        for i=1:nz
            for j=1:ne
                Err=Err+abs(sfold(i,j)-sf(i,j));
            end;
        end;
        if Err <= MaxErr
            break;
        end;%stop if iteration has converged
    end;
    Err,iter
    %define boundary conditions for vorticity
    for i=1:nz
        vt(i,1)=-2*(xz(i,1)^2+yz(i,1)^2)*(sf(i,1)-sf(i,2))/J(i,1)^2;
        vt(i,ne)=-2*(xz(i,ne)^2+yz(i,ne)^2)*(sf(i,ne)-sf(i,ne-1))/J(i,ne)^2;
    end;
    for j=2:ne-1
        vt(1,j)=(6*Q/((ne-1)*h)^2)*(1-2*((j-1)/(ne-1)));
        vt(nz,j)=vt(nz-1,j);
    end;
    %calculate right hand side of vorticity
    for i=2:nz-1
        for j=2:ne-1
            R(i,j) = q1(i,j)*(vt(i+1,j)-2*vt(i,j)+vt(i-1,j))/h^2/J(i,j)^2-...
                     2*q2(i,j)*(vt(i+1,j+1)-vt(i-1,j+1)-vt(i+1,j-1)+vt(i-1,j-1))/4/h^2/J(i,j)^2+...
                     q3(i,j)*(vt(i,j+1)-2*vt(i,j)+vt(i,j-1))/h^2/J(i,j)^2+...
                     Lz(i,j)*(vt(i+1,j)-vt(i-1,j))/2/h+...
                     Le(i,j)*(vt(i,j+1)-vt(i,j-1))/2/h;
            w(i,j) = alpha*R(i,j)-...
                     (sf(i,j+1)-sf(i,j-1))*(vt(i+1,j)-vt(i-1,j))/4/h^2/J(i,j)+...
                     (sf(i+1,j)-sf(i-1,j))*(vt(i,j+1)-vt(i,j-1))/4/h^2/J(i,j);
        end;
    end;
    %update vorticity
    vt(2:nz-1,2:ne-1)=vt(2:nz-1,2:ne-1)+dt*w(2:nz-1,2:ne-1);
    
    %update time
    t=t+dt;
%plot
%figure(2);
pause(0.1);
hold on;
subplot(211);contourf(x,y,vt);title('Vorticity');
subplot(212);contourf(x,y,sf);title('Stream Function');
hold off;
end;

        








