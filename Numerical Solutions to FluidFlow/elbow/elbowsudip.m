clear;clc;
gridgen;
dxi = 1.0; 
deta = 1.0; 
q = 1.0; alpha = 0.1; t = 0; dt = 0.01; max_step = 100; max_iter = 1000; 
tol = 1.0e-03; beta = 1.5;

disp('SOLVE FOR STREAM FUNCTION INSIDE DOMAIN') 

sf = zeros(NTot,MTot); sf_old = zeros(NTot,MTot); vt = zeros(NTot,MTot); rhs = zeros(NTot,MTot); 
sf(1:NTot,1) = 0.0;%sf at lower wall 
sf(1:NTot,MTot) = q;%sf at upper wall 
for j = 2:MTot-1; sf(1,j) = q*(3*((j-1)/(MTot-1))^2-2*((j-1)/(MTot-1))^3); end; %sf at inlet 
disp('START TIME INTEGRATION') 
for tstep = 1:max_step,tstep,t 
error = 999; iter = 0; 
sf
while error > tol 
sf(NTot,2:MTot-1) = sf(NTot-1,2:MTot-1);%sf at outlet, derivative = 0 
sf_old = sf; error = 0; iter = iter + 1;
for i = 2:NTot-1; for j = 2:MTot-1 
dxdxi = (x(i+1,j)-x(i-1,j))/(2*dxi); 
dxdeta = (x(i,j+1)-x(i,j-1))/(2*deta); 
d2xdxi2 = (x(i+1,j)-2*x(i,j)+x(i-1,j))/(dxi^2); 
d2xdeta2 = (x(i,j+1)-2*x(i,j)+x(i,j-1))/(deta^2); 
d2xdxideta = (x(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1))/(4*dxi*deta); 
dydxi = (y(i+1,j)-y(i-1,j))/(2*dxi); 
dydeta = (y(i,j+1)-y(i,j-1))/(2*deta); 
d2ydxi2 = (y(i+1,j)-2*y(i,j)+y(i-1,j))/(dxi^2); 
d2ydeta2 = (y(i,j+1)-2*y(i,j)+y(i,j-1))/(deta^2); 
d2ydxideta = (y(i+1,j+1)-y(i+1,j-1)-y(i-1,j+1)+y(i-1,j-1))/(4*dxi*deta); 
J = dxdxi*dydeta-dxdeta*dydxi; 
q1 = dxdeta^2+dydeta^2; 
q2 = dxdxi*dxdeta+dydxi*dydeta; 
q3 = dxdxi^2+dydxi^2; 
del2xi = (1/J^3)*(q1*(dxdeta*d2ydxi2-dydeta*d2xdxi2)-2*q2*(dxdeta*d2ydxideta-... 
dydeta*d2xdxideta)+q3*(dxdeta*d2ydeta2-dydeta*d2xdeta2)); 
del2eta = (1/J^3)*(q1*(dydxi*d2xdxi2-dxdxi*d2ydxi2)-2*q2*(dydxi*d2xdxideta-... 
dxdxi*d2ydxideta)+q3*(dydxi*d2xdeta2-dxdxi*d2ydeta2)); 
a0 = (2/J^2)*(q1/dxi^2+q3/deta^2); 
a1 = q1/(J^2*dxi^2)+del2xi/(2*dxi); 
a2 = q1/(J^2*dxi^2)-del2xi/(2*dxi); 
a3 = q3/(J^2*deta^2)+del2eta/(2*deta); 
a4 = q3/(J^2*deta^2)-del2eta/(2*deta); 
a5 = q2/(2*J^2*dxi*deta); 
sf(i,j) = beta*(1/a0)*( (a1*sf(i+1,j)+a2*sf(i-1,j)+a3*sf(i,j+1)+a4*sf(i,j-1))-... 
a5*(sf(i+1,j+1)-sf(i+1,j-1)-sf(i-1,j+1)+sf(i-1,j-1))+... 
vt(i,j) ) + (1-beta)*sf_old(i,j); 
error = error + (abs(sf(i,j))-abs(sf_old(i,j)))^2; 
end; end; 
error = sqrt(error); 
if iter >= max_iter, break; end 
end 
sf(NTot,2:MTot-1) = sf(NTot-1,2:MTot-1);%sf at outlet, derivative = 0 
iter 

disp('UPDATE BOUNDARY VORTICITY') 

for i = 2:NTot-1
j=MTot-1;
dxdxi = (x(i+1,j)-x(i-1,j))/(2*dxi); 
dxdeta = (x(i,j+1)-x(i,j-1))/(2*deta); 
dydxi = (y(i+1,j)-y(i-1,j))/(2*dxi); 
dydeta = (y(i,j+1)-y(i,j-1))/(2*deta); 
J = dxdxi*dydeta-dxdeta*dydxi; 
vt(i,1) = (2*(dxdxi^2+dydxi^2)/J^2)*(sf(i,1)-sf(i,2));%vt at lower wall 
vt(i,MTot) = (-2*(dxdxi^2+dydxi^2)/J^2)*(sf(i,MTot)-sf(i,MTot-1));%vt at upper wall 
end 
for i = 2:NTot-1
j=2;
dxdxi = (x(i+1,j)-x(i-1,j))/(2*dxi); 
dxdeta = (x(i,j+1)-x(i,j-1))/(2*deta); 
dydxi = (y(i+1,j)-y(i-1,j))/(2*dxi); 
dydeta = (y(i,j+1)-y(i,j-1))/(2*deta); 
J = dxdxi*dydeta-dxdeta*dydxi; 
vt(i,1) = (2*(dxdxi^2+dydxi^2)/J^2)*(sf(i,1)-sf(i,2));%vt at lower wall 
vt(i,MTot) = (-2*(dxdxi^2+dydxi^2)/J^2)*(sf(i,MTot)-sf(i,MTot-1));%vt at upper wall 
end 
for j = 2:MTot-1; vt(1,j) = (6*q/0.5^2)*(1-2*j/MTot); end;%vt at inlet 
vt(NTot,2:MTot-1) = vt(NTot-1,2:MTot-1);%vt at outlet, derivative = 0 

disp('COMPUTE RHS OF VORTICITY EQUATION') 
for i = 2:NTot-1; for j = 2:MTot-1 
dxdxi = (x(i+1,j)-x(i-1,j))/(2*dxi); 
dxdeta = (x(i,j+1)-x(i,j-1))/(2*deta); 
d2xdxi2 = (x(i+1,j)-2*x(i,j)+x(i-1,j))/(dxi^2); 
d2xdeta2 = (x(i,j+1)-2*x(i,j)+x(i,j-1))/(deta^2); 
d2xdxideta = (x(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1))/(4*dxi*deta); 
dydxi = (y(i+1,j)-y(i-1,j))/(2*dxi); 
dydeta = (y(i,j+1)-y(i,j-1))/(2*deta); 
d2ydxi2 = (y(i+1,j)-2*y(i,j)+y(i-1,j))/(dxi^2); 
d2ydeta2 = (y(i,j+1)-2*y(i,j)+y(i,j-1))/(deta^2); 
d2ydxideta = (y(i+1,j+1)-y(i+1,j-1)-y(i-1,j+1)+y(i-1,j-1))/(4*dxi*deta); 
J = dxdxi*dydeta-dxdeta*dydxi; 
q1 = dxdeta^2+dydeta^2; 
q2 = dxdxi*dxdeta+dydxi*dydeta; 
q3 = dxdxi^2+dydxi^2; 
del2xi = (1/J^3)*(q1*(dxdeta*d2ydxi2-dydeta*d2xdxi2)-2*q2*(dxdeta*d2ydxideta-... 
dydeta*d2xdxideta)+q3*(dxdeta*d2ydeta2-dydeta*d2xdeta2)); 
del2eta = (1/J^3)*(q1*(dydxi*d2xdxi2-dxdxi*d2ydxi2)-2*q2*(dydxi*d2xdxideta-... 
dxdxi*d2ydxideta)+q3*(dydxi*d2xdeta2-dxdxi*d2ydeta2)); 
sfxi = (sf(i+1,j)-sf(i-1,j))/(2*dxi); 
sfeta = (sf(i,j+1)-sf(i,j-1))/(2*deta); 
vtxi = (vt(i+1,j)-vt(i-1,j))/(2*dxi); 
vteta = (vt(i,j+1)-vt(i,j-1))/(2*deta); 
vt2dxi2 = (vt(i+1,j)-2*vt(i,j)+vt(i-1,j))/(dxi^2); 
vt2deta2 = (vt(i,j+1)-2*vt(i,j)+vt(i,j-1))/(deta^2); 
vt2dxideta = (vt(i+1,j+1)-vt(i+1,j-1)-vt(i-1,j+1)+vt(i-1,j-1))/(4*dxi*deta); 
rhs(i,j) = -(1/J)*(sfeta*vtxi-sfxi*vteta)+alpha*((1/J^2)*(q1*vt2dxi2-... 
2*q2*vt2dxideta+q3*vt2deta2)+del2xi*vtxi+del2eta*vteta); 
end; end; 

disp('SOLVE FOR VORTICITY INSIDE DOMAIN') 

vt(2:NTot-1,2:MTot-1) = vt(2:NTot-1,2:MTot-1)+dt*rhs(2:NTot-1,2:MTot-1); 

t = t + dt;%next time step 

hold on;
subplot(211);contourf(x,y,vt);title('Vorticity');
subplot(212);contourf(x,y,sf);title('Stream Function');
hold off;
pause(0.1);

end 
