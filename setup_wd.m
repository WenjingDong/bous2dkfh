clear all;
close all;
nx = 256;
nz = 128;
L = 2*pi; 
H = pi;
dx = L/nx; 
dz = H/nz;
x = (0:nx-1)*dx;  
z = (0:nz)*dz; % functions are not periodic in z
[x_,z_] = ndgrid(x,z); 

N = 15;
f = 0.;
k = 6;
m = 3;
kappa = sqrt(k^2+m^2);
omega = sqrt(k^2*N^2 + m^2*f^2)/kappa;
a = .05;
factor = 1;

%psi = a*cos(k*x_).*sin(m*z_);
% If I added the second term, u blows up
w =  a*k*sin(k*x_).*sin(m*z_);
V =   a*f^2*m/omega * sin(k*x_).*cos(m*z_);
b =  -N^2*a*k/omega * cos(k*x_).*sin(m*z_);

u =  a*m*cos(k*x_).*cos(m*z_) + factor * k/(2*omega)*(a*k)^2 * (m^2/k^2+1) * (1 -  cos(2*m*z_))/2;




figure;
pcolor(x,z,u');shading interp;colorbar;

Sin(:,:,1) = u;
Sin(:,:,2) = w;
Sin(:,:,3) = V;
Sin(:,:,4) = b;

% increasing the number of particles is not very helpful
np = 64;
x0 = 2*pi/np/2;
Xpin.y = pi/2*ones(np,1); 
Xpin.x = linspace(x0,2*pi+x0-2*pi/np,np)';

%figure
%pcolor(x,z,b'), shading interp, colorbar, axis image

numsteps = 50000;
savestep = 100;
[Sout,time,energy,Xp] = bous2kfh(Sin,N,f,numsteps,savestep,Xpin);

% validation of the code: compare the velocity with exact soln
u_ex = a*m*cos(k*x_ - omega * time(end)).*cos(m*z_);
w_ex = a*k*sin(k*x_ - omega * time(end)).*sin(m*z_);
b_ex = -N^2*a*k/omega * cos(k*x_ - omega * time(end)).*sin(m*z_);

%figure;
%contourf(x,z,w_ex' - squeeze( Sout(:,:,2,end))'); colorbar;

xpar = [Xp(:).x];
xpar_disp = xpar - xpar(:,1);
zpar = [Xp(:).y];
zpar_disp = zpar- zpar(:,1);
for j=1:numsteps+1
    Xpa(j) = mean(xpar_disp(:,j));
    zpa(j) = mean(zpar_disp(:,j));
end



%us = a^2*k*m^2/(2*omega);
us = k/(2*omega) * m^2/k^2 * (a*k)^2 *cos(2*m*pi/2);
ul = k/(2*omega) * (a*k)^2 * (m^2/k^2 * (cos(m*pi/2))^2 + (sin(m*pi/2))^2);



dt = (time(end)-time(end-1))/savestep;
t=[0:numsteps]*dt;

figure('units','centimeter','Position', [10 5 10 8]);
%subplot(1,2,1);
plot(t,Xpa,'b','linewidth',2);hold on;
plot(t,us*t,'r','linewidth',2);
plot(t,ul*t,'g','linewidth',2);

%plot(t(2:end),xpar_disp(2,2:end),'k');
%plot(t(2:end),xpar_disp(5,2:end),'m');

lg = legend('$<x>$', '$\overline{u}^st$','$\overline{u}^Lt$');
set(lg,'interpreter','latex','FontSize',12);
legend boxoff;
xlabel('t'); ylabel('x');
set(gca,'FontSize',12);
%title('particle displacement in x');
grid on;
title('m=2')
subplot(1,2,2);
plot(t(2:end),(Xpa(2:end) - ul*t(2:end))./(ul*t(2:end)) ,'b','linewidth',2);
hold on;
plot(t(2:end),(us*t(2:end) - ul*t(2:end))./(ul*t(2:end)),'r','linewidth',2);
%plot(t,ul*t,'r','linewidth',2);

lg = legend('$\frac{<x>-\overline{u}^Lt}{\overline{u}^Lt}$', '$\frac{\overline{u}^st-\overline{u}^Lt}{\overline{u}^Lt}$');
set(lg,'interpreter','latex','FontSize',14);
legend boxoff;
xlabel('t'); 
set(gca,'FontSize',12);
grid on;



print -f3 -depsc -r600 k6m2_f0

%plot(t,xpar_disp(9,:)- ul*t,'k--');
%plot(t,xpar_disp(10,:)- ul*t,'b');
 %plot(xlim,[0 0],'k--');

% yay!
