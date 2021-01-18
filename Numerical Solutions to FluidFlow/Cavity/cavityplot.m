%Plotting Velocity and Vorticity
figure(1);hold off;
quiver(flipud(rot90(uu)),flipud(rot90(vv)),'r');
hold on;
contour(flipud(rot90(w)),20);axis equal;
title(sprintf('Steady State Velocity and Vorticity Plot for %d X %d grid',nx,ny));
xlabel('X');ylabel('Y');
%Plotting Pressure
figure(2);hold off;
x=1:nx;y=1:ny;
surf(x,y,pp);axis equal;
title(sprintf('Steady State Pressure Plot for %d X %d grid',nx,ny));
xlabel('X');ylabel('Y');zlabel('Pressure');
%Plotting Stream Function
figure(3);hold off;
contour(flipud(rot90(sf)),20);axis equal;
pause(0.001);
title(sprintf('Steady State Stream function Plot for %d X %d grid',nx,ny));
xlabel('X');ylabel('Y');
