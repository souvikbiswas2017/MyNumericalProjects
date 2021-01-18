%-----------------------------------------------------------
%Start Grid Generation
%-----------------------------------------------------------
nB=2;%Number of blocks
%lower boundary
xp(1,1:nB+1)=[0,0,1.5];yp(1,1:nB+1)=[1,0,0];
%upper boundary
xp(2,1:nB+1)=[0.5,0.5,1.5];yp(2,1:nB+1)=[1,0.5,0.5];

Nblock(1:nB)=[20,20];NTot=0;
Mblock(1:nB)=[15,15];MTot=15;

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
for m=1:1
	for i=2:NTot-1
	for j=2:MTot-1
		x(i,j)=0.25*(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1));
		y(i,j)=0.25*(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1));
	end;
	end;
end;
title('Grid for elbow');
axis('square'),axis([-0.10, 2, 0, 2]); hold on
for j=1:MTot;plot(x(1:NTot,j),y(1:NTot,j),'r','linewidth',2);end;
for i=1:NTot;plot(x(i,1:MTot),y(i,1:MTot),'r','linewidth',2);end;
axis off
%-----------------------------------------------------------
%End Grid Generation
%-----------------------------------------------------------
