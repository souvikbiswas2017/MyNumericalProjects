clf;nx=17; ny=17;MaxStep=120; Visc=0.1; dt=0.01;         %resolution & parameters governing the run
MaxIt=100; Beta=1.5; MaxErr=0.001;                    %parameters for SOR iteration
sf=zeros(nx,ny); vt=zeros(nx,ny); w=zeros(nx,ny); h=1.0/(nx-1); t=0.0;
for istep=1:MaxStep,                                  %start the time integration
  for iter=1:MaxIt,                                   %solve for the streamfunction
    w=sf;                                             %by SOR iteration
    for i=2:nx-1; for j=2:ny-1
        sf(i,j)=0.25*Beta*(sf(i+1,j)+sf(i-1,j)...
          +sf(i,j+1)+sf(i,j-1)+h*h*vt(i,j))+(1.0-Beta)*sf(i,j);
       end; end;
    Err=0.0;
    for i=1:nx,for j=1:ny, Err=Err+abs(w(i,j)-sf(i,j)); end; end;
    if Err <= MaxErr, break, end                      %stop if iteration has converged
  end;
    vt(2:nx-1,1)=-2.0*sf(2:nx-1,2)/(h*h);                    %set vorticity on bottom wall
    vt(2:nx-1,ny)=-2.0*sf(2:nx-1,ny-1)/(h*h)-2.0/h;                %top wall
    vt(1,2:ny-1)=-2.0*sf(2,2:ny-1)/(h*h);                          %right wall
    vt(nx,2:ny-1)=-2.0*sf(nx-1,2:ny-1)/(h*h);                      %left wall
  for i=2:nx-1; for j=2:ny-1                                       % compute the
      w(i,j)=-0.25*((sf(i,j+1)-sf(i,j-1))*(vt(i+1,j)-vt(i-1,j))... % right hand side of the 
       -(sf(i+1,j)-sf(i-1,j))*(vt(i,j+1)-vt(i,j-1)))/(h*h)...      % vorticity equation
       +Visc*(vt(i+1,j)+vt(i-1,j)+vt(i,j+1)+vt(i,j-1)-4.0*vt(i,j))/(h*h);
    end; end;
  vt(2:nx-1,2:ny-1)=vt(2:nx-1,2:ny-1)+dt*w(2:nx-1,2:ny-1);         %update the vorticity
  t=t+dt                                                           %print out t every step
%   subplot(121), contour(rot90(fliplr(vt))), axis('square');        %plot vorticity
%   subplot(122), contour(rot90(fliplr(sf))), axis('square');pause(0.01) %streamfunction
  subplot(121), contourf(rot90(fliplr(vt))), axis('square');        %plot vorticity
  subplot(122), contourf(rot90(fliplr(sf))), axis('square');pause(0.01) %streamfunction
end;
