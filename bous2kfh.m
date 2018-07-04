function [Sout,time,energy,Xp,bmov,Bp] = bous2kfh(Sin,N,f,numsteps,savestep,Xpin)

%  [Sout,time,energy,Xp,bmov] = bous2kfm(Sin,N,f,numsteps,savestep,Xpin)
%
%  Solves 2D Boussinesq equations in x-z plane 
%
%  Du/Dt - V + P_x = d(u)
%  Dw/Dt - b + P_z = d(w)
%  Dv/Dt + f u = d(v)       use v instead of V=fv to avoid divide V by zero
%  Db/Dt + N^2 w = d(b)
%  del^2 P = V_x + b_z + 2J(u,w)
%
%  J(u,w) =  u_x w_z - u_z w_x  
%  D/Dt = \p_t + u\p_x + w\p_z
%  d(A) =  nu*del^(hv) A
%  (u,w) = (psi_z,-psi_x)
%  zeta = psi_xx + psi_zz
%
%  Inputs
%
%  Sin:       nx x nz+1 x 4 array with initial u, w, V, b (resp)
%  The dimension of Sin has been increased to include the top level, z=H. (added by Wenjing)
%  N:         Mean buoyancy gradient (H/pi * N^2_dim)
%  f:         Coriolis parameter
%  numsteps:  Total number of timesteps 
%  savestep:  Frequency, in timesteps, to save output    
%  np:        Advect np^2 particles (optional)
%
%  Outputs
%
%  Sout:      Arranged as Sin, but with 4th dimension for time
%  time:      Times at which output is saved
%  energy:    Time series of energy
%  bmov:      Movie of b field (if bmov included in output list)
%  Xp:        Structure with coordinates Xp.x, Xp.y, Xp.z of particles
%    
%  Numerical details
%
%  Model is spectral, in square domain of size 2*pi x pi.  Input fields
%  must have nx = 2^n, where n is an integer, and likewise for nz.
%  Nonlinear terms are done in physical space using dealiased product
%  in x, but not z.  Uses sine and cosine transforms, as appropriate,
%  for z.  Uses AB3 timestepping with trapezoidal hyperviscosity of
%  order a.  Timestep and hyperviscosity are set adaptively, via dt =
%  dttune*dx/max(|u|) and nu = nutune*dx^hv*zeta_rms.
%
%  sin and cos transforms done using fft, with even or odd
%  symmmetry.  u, V are even;  w, b are odd;  ddz <-> ikz_ switches symmetry

%  ------------------
%   Added by Wenjing
%  When a field variable is even, its value at z=H has to be specified rather than 
%  calculated by forcing the average to be zero.
%  Thus the dimension of the data is nx * (nz+1)
%  ------------------

%  Tuning factors dttune and nutune and hyperviscous order a
%  can be set by editing this file directly.

% Tuning
hv = 8;        % hyperviscosity order, nu del^hv u 
nutune = 10;
dttune = .2;   % Courant number
   
% Set global params for use in rhs functions
global nx nz ikx_ ikz_ K2_ da noda odd even u_ind w_ind V_ind b_ind

N2 = N^2;  %f2 = f^2;

u_ind = 1; w_ind = 2; v_ind = 3; b_ind = 4;

if nargin > 5, 
    np = length(Xpin.x);
    Xp = Xpin;
else 
    np = 0;
    Xp = 0;
end

% Check for outputs requested
makemov = false;
if (nargout>4), makemov=true;  end

% Get and check dimensions
[nx,nzt,nS] = size(Sin);
% nz' is a temporary variable
nz = nzt - 1;

if (mod(log2(nx),1)~=0), error('must have nx = 2^(integer)'); end
if (mod(log2(nz),1)~=0), error('must have nz = 2^(integer)'); end

% Set grid
dx = 2*pi/nx;
dz = pi/nz;
x = (0:nx-1)*dx;  
z = (0:nz)*dz;
[x_,z_] = ndgrid(x,z); 

% Set up arrays of wavenumbers for computing derivatives
[ikx_,ikz_] = ndgrid(0:nx/2-1,0:nz-1); % nx/2 x nz
K2_   = ikx_.^2 + ikz_.^2;
K2_(1,1) = .1;  % prevent divide-by-zero in computation of psi
ikx_ = sqrt(-1)*ikx_; 
ikz_ = sqrt(-1)*ikz_;

% Flags for transforms
da  = 1;  noda = 0;
odd = 1;  even = 0;

% Trapezoidal hyperdiffusion operators
nudt = nutune*2*pi/(nx*(nx/2-1)^hv);  % divide by dt to get effective nu
fR = (1+nudt/2*sqrt(K2_).^hv).^(-1);
fU = (1-nudt/2*sqrt(K2_).^hv).*fR;
filterU = repmat(fU,[1 1 nS]);       
filterR = repmat(fR,[1 1 nS]);
clear fR fU

% Get initial spectral PV - both b and zeta have odd-symmetry in z

size(Sin)
Sk(:,:,u_ind) = g2k(mirror(Sin(:,:,u_ind),even));  % u is even
Sk(:,:,w_ind) = g2k(mirror(Sin(:,:,w_ind),odd));   % w is odd
Sk(:,:,v_ind) = g2k(mirror(Sin(:,:,v_ind),even));  % V is even
Sk(:,:,b_ind) = g2k(mirror(Sin(:,:,b_ind),odd));   % b is odd

figure;
pcolor(mirror(Sin(:,:,u_ind),even)');shading interp;colorbar;
% Maximum phase speed.  Note c_x = +/- N/sqrt(k^2+m^2)
% gravest modes:  k = 1, m = 1/2
Cmax = sqrt(N2);

% Set params for AB3+trapezoidal diffusion (Durran 3.81)
a1 = 23/12;   a2 = -16/12;  a3 = 5/12;  

% Preallocate arrays to hold saved output
nframes = floor(numsteps/savestep);
Sout = zeros(nx,nz+1,nS,nframes);  % output on grid
time = zeros(1,nframes); 
ke = zeros(1,nframes); 
pe = zeros(1,nframes); 
if (np>0)
    Bc = interpolate(Xpin.x,Xpin.z,Sin(:,:,b_ind)+N2*z_,dx,dz);
    Bc = [0; Bc];
end

% Set counters and temp fields
frame = 0;  t = 0;  n = 0; 
Rk = 0; Rkm1 = 0; Rkm2 = 0;

keepgoing = true;
while keepgoing
        
    % Save n-1 and n-2 rhs and get next one.  
    Rkm2 = Rkm1;
    Rkm1 = Rk;
    Rk = getrhs(Sk,N2,f);  
    
    if (n==0) Rkm1 = Rk; Rkm2 = Rk; end

    u = k2g(Sk(:,:,u_ind),even,noda);
    w = k2g(Sk(:,:,w_ind),odd,noda);
    Umax = max([max(abs(real(u(:)))) max(abs(real(w(:)))) Cmax]);
    
    % Exit if blowing up, and save last field.
    if (Umax>1e6)  
        disp(strcat('Blow up! Umax= ',num2str(Umax),', t= ',num2str(t)))
        Sout(:,:,u_ind,frame+1) = k2g(Sk(:,:,u_ind),even,noda); 
        Sout(:,:,w_ind,frame+1) = k2g(Sk(:,:,w_ind),odd,noda); 
        Sout(:,:,v_ind,frame+1) = k2g(Sk(:,:,v_ind),even,noda); 
        Sout(:,:,b_ind,frame+1) = k2g(Sk(:,:,b_ind),odd,noda); 
        keepgoing = false;
    end
    
    % Adapt dt and nu
    dt = dttune*dx/Umax;   % Courant condition 
    nu = nudt/dt;
        
    % Save output at frequency savestep
    if (mod(n,savestep)==0||n==0)  
        frame = frame+1;  
        u = k2g(Sk(:,:,u_ind),even,noda); 
        w = k2g(Sk(:,:,w_ind),odd,noda); 
        v = k2g(Sk(:,:,v_ind),even,noda); 
        b = k2g(Sk(:,:,b_ind),odd,noda); 
        
        % check if the flow is non-divergent
        div = k2g(ikx_.*Sk(:,:,u_ind)+ikz_.*Sk(:,:,w_ind),even,noda); 
        
        Sout(:,:,u_ind,frame) = u;
        Sout(:,:,w_ind,frame) = w;
        Sout(:,:,v_ind,frame) = v;
        Sout(:,:,b_ind,frame) = b;
        
        divtotal(frame) = sum(abs(div(:)));
        energy(frame) = .5*sum(u(:).^2+(V(:)/f).^2+w(:).^2+b(:).^2/N2)/(nx*nz);
        time(frame) = t;
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        disp(strcat('max(|u|) = ',num2str(Umax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        disp(strcat('divergent: ', num2str(divtotal(frame))));
        if (makemov)
            hc = figure('units','centimeters','position',[5 5 10 30]);
            clf
            subplot(3,1,1)
            %contourf(x,z,(b+N2*z_)',Bc), colorbar, axis image; 
            contourf(x,z,(b+N2*z_)'), colorbar, axis image; 
            title('b + N^2z')
            xlabel('x')
            ylabel('z')
            set(gca,'fontsize',16)
            if (np>0)
                hold
                xp=[Xp(n+1).x]; yp=[Xp(n+1).y];
                plot(xp,yp,'k.','MarkerSize',10)
            end
            subplot(3,1,2)
            %psi = k2g(-Sk(:,:,1)./K2_,odd,noda); 
            
            % added by Wenjing. calcualte PV -fN^2
            b_z = k2g(ikz_.*Sk(:,:,b_ind),even,noda);
            v_x =  k2g(ikx_.*Sk(:,:,v_ind),even,noda); 
            b_x = k2g(ikx_.*Sk(:,:,b_ind),odd,noda);
            v_z =  k2g(ikz_.*Sk(:,:,v_ind),odd,noda);
            
            u_z = k2g(ikz_.*Sk(:,:,u_ind),odd,noda); %% this one has a problem
            w_x = k2g(ikx_.*Sk(:,:,w_ind),odd,noda);
            
            q = f*b_z + N^2 * v_x + b_x .* (-v_z) + b_z .* v_x;
            zeta = u_z - w_x;
            
            pcolor(x,z,q'), shading interp, colorbar, axis image;
            title('Vorticity q')
            set(gca,'fontsize',16)
            

            max(max(Sk(:,:,u_ind)))
            subplot(3,1,3);
            pcolor(x,z,u');shading interp, colorbar, axis image;
            drawnow
            bmov(frame) = getframe(hc);
        end
    end
    
    % Timestep and diffuse
    Sk = filterU.*Sk + dt*filterR.*(a1*Rk + a2*Rkm1 + a3*Rkm2);
    
    %%% added by Wenjing
 
    
    

    
    if (np>0) % advect particles -- MAKE 3D?!
        Bp(:,n+1) = interpolate(Xp(n+1).x,Xp(n+1).z,b+N2*z_,dx,dz);
        Xp(n+2) = advect_particles(Xp(n+1),u,v, w,dx,dz,dt);
        % Get total buoyancy b+N^2*z at particle locations
    end
    
    n = n+1;
    t = t+dt;  % clock
    
    if (n==numsteps), disp('End reached'), keepgoing=false; end

     
end

return

%-------------------------------------------------------------------
% Internal functions
%-------------------------------------------------------------------

function Rk = getrhs(Sk,N2,f)
    
    %  u_t = -u u_x - w u_z + V - P_x  (u even)
    %  w_t = -u w_x - w w_z + b - P_z  (w odd)
    %  v_t = -u v_x - w v_z - f u    (V even)
    %  b_t = -u b_x - w b_z - N^2 w    (b odd)
    %  del^2 P = 2[u_x w_z - u_z w_x] + fv_x + b_z 

    global ikx_ ikz_ K2_ odd even da u_ind w_ind V_ind b_ind
    
    u   = k2g(Sk(:,:,u_ind),even,da);
    w   = k2g(Sk(:,:,w_ind),odd,da);
    u_x = k2g(ikx_.*Sk(:,:,u_ind),even,da);
    w_x = k2g(ikx_.*Sk(:,:,w_ind),odd,da);
    v_x = k2g(ikx_.*Sk(:,:,v_ind),even,da);
    b_x = k2g(ikx_.*Sk(:,:,b_ind),odd,da);
    u_z = k2g(ikz_.*Sk(:,:,u_ind),odd,da);
    w_z = k2g(ikz_.*Sk(:,:,w_ind),even,da);
    v_z = k2g(ikz_.*Sk(:,:,v_ind),odd,da);
    b_z = k2g(ikz_.*Sk(:,:,b_ind),even,da);
    
    Pk = -(2*g2k(u_x.*w_z - u_z.*w_x) + ikx_.*f*Sk(:,:,v_ind) + ikz_.*Sk(:,:,b_ind))./K2_;
    Rk(:,:,u_ind) = -g2k(u.*u_x + w.*u_z) + Sk(:,:,V_ind) - ikx_.*Pk;
    Rk(:,:,w_ind) = -g2k(u.*w_x + w.*w_z) + Sk(:,:,b_ind) - ikz_.*Pk;  
    Rk(:,:,v_ind) = -g2k(u.*v_x + w.*v_z) - f*Sk(:,:,u_ind);  
    Rk(:,:,b_ind) = -g2k(u.*b_x + w.*b_z) - N2*Sk(:,:,w_ind);  
    
    return
    
%-------------------------------------------------------------------
% $$$ 
% $$$ function Pk = getP(Sk)
% $$$ 
% $$$ global ikx_ ikz_ K2_ odd even da u_ind w_ind V_ind b_ind
% $$$ 
% $$$ u_x = k2g(ikx_.*Sk(:,:,u_ind),even,da);
% $$$ w_x = k2g(ikx_.*Sk(:,:,w_ind),odd,da);
% $$$ u_z = k2g(ikz_.*Sk(:,:,u_ind),odd,da);
% $$$ w_z = k2g(ikz_.*Sk(:,:,w_ind),even,da);
% $$$ 
% $$$ Pk = -(2*g2k(u_x.*w_z - u_z.*w_x) + ikx_.*Sk(:,:,V_ind) + ikz_.*Sk(:,:,b_ind))./K2_;

%-------------------------------------------------------------------

function fgf = mirror(fg,sym)

    % Populate upper-half-plane in z with odd or even mirrored field
    % Input:  fg(nx,nz)
    % Output: fgm(nx,2*nz)
    % sym = 1 => odd, sym = 0 => even.

    %---------
    % values at z=H are given by initial conditions
    %---------       
    global nx nz
    
    fgf = zeros(nx,2*nz);
    fgf(:,1:nz+1) = fg;

    if sym % odd
%        fgf(:,nz+1) = fg(:,1);
        fgf(:,end:-1:nz+2) = -fg(:,2:end-1);
    else % even
        fgf(:,end:-1:nz+2) = fg(:,2:end-1);
%        fgf(:,nz+1) = -sum(fgf,2); % demand mean 0 .. nx/2+1 still 0
    end
    
    return
        
%-------------------------------------------------------------------

function fk = g2k(fg)
    
    % Input:  fg(n1,n2)
    % Output: fk(nx/2,nz)
    
    global nx nz
    
    % n1 should be either 3*nx/2 or nx.  n2 should be 2*nz.
    [n1,n2] = size(fg);
    fkf = fft2(fg)/(n1*n2);
    fk = fkf(1:nx/2,1:nz);
    
    return
        
%-------------------------------------------------------------------

function fg = k2g(fk,sym,da)

    % Input:  fk(nx/2,nz)
    % Transform to grid space, demanding conjugate symmetry 
    % and odd or even symmetry in second coordinate.  
    % To use for nonlinear products,
    % set da = 1 => Output:  fg(3*nx/2,2*nz).  
    % For saving output as gridded fields, 
    % set da = 0 => Output:  fg(nx,2*nz) 
    % sym = 1 => odd, sym = 0 => even.

    global nx nz
    
    fkf = zeros(nx,2*nz);
    fkf(1:nx/2,1:nz) = fk;
    if sym % odd
        fkf(end:-1:nx/2+2,1:nz) = -conj(fkf(2:nx/2,1:nz));
        fkf(:,end:-1:nz+2) = -fkf(:,2:nz);
    else % even
        fkf(end:-1:nx/2+2,1:nz) = conj(fkf(2:nx/2,1:nz));
        fkf(:,end:-1:nz+2) = fkf(:,2:nz);
    end
    
    if da, nxh = 3*nx/2; else nxh = nx; end
    
    fg = nxh*2*nz*ifft(ifft(fkf,2*nz,2),nxh,'symmetric'); % nxh x 2*nz
   
    % keep values at level nz+1, i.e. z=H. 
    if ~da, fg = fg(:,1:nz+1); end
        
    fg = real(fg);
    
    return
        
%-------------------------------------------------------------------

function [Xp] = advect_particles(Xp,u,v,w, dx,dz,dt)
    
% particle position vector for particle m: Xp(t)=[Xp(t).x_m,Xp(t).y_m]
% m = 1:np
% Inputs:  Xp:  structure with Xp.x and Xp.y
%          u, v, w:  velocity field
    
% advect with RK4

xp1 = dt*interpolate(Xp.x,Xp.z,u,dx,dz);
yp1 = dt*interpolate(Xp.x,Xp.z,v,dx,dz);
zp1 = dt*interpolate(Xp.x,Xp.z,w,dx,dz);

xp2 = dt*interpolate(Xp.x+xp1/2,Xp.z+zp1/2,u,dx,dz);
yp2 = dt*interpolate(Xp.x+xp1/2,Xp.z+zp1/2,v,dx,dz);
zp2 = dt*interpolate(Xp.x+xp1/2,Xp.z+zp1/2,w,dx,dz);

xp3 = dt*interpolate(Xp.x+xp2/2,Xp.z+zp2/2,u,dx,dz);
yp3 = dt*interpolate(Xp.x+xp2/2,Xp.z+zp2/2,v,dx,dz);
zp3 = dt*interpolate(Xp.x+xp2/2,Xp.z+zp2/2,w,dx,dz);

xp4 = dt*interpolate(Xp.x+xp3,Xp.z+zp3,u,dx,dz);
yp4 = dt*interpolate(Xp.x+xp3,Xp.z+zp3,v,dx,dz);
zp4 = dt*interpolate(Xp.x+xp3,Xp.z+zp3,w,dx,dz);

Xp.x = Xp.x + (xp1 + 2*xp2 + 2*xp3 + xp4)/6;
Xp.y = Xp.y + (yp1 + 2*yp2 + 2*yp3 + yp4)/6; 
Xp.z = Xp.z + (zp1 + 2*zp2 + 2*zp3 + zp4)/6; 

%-------------------------------------------------------------------

function FI = interpolate(x, y, F, dx, dy)

% Interpolate function F (velocity component u or v, on grid) to
% particle positions (x,y). Uses Lagrangian interpolation of order
% Iord (Iord = 1 ==> cubic interpolation).  See
% Duran Ch. 6, for example.  ax,ay below are x and y components
% of what he calls alpha, the fractional grid position.

% nx = 2*(kmax+1)
% dx = 2*pi/nx 

Iord = 1;
bump = 10^(-10); % Prevent NaNs

[nx,ny] = size(F);
FI = zeros(size(x));  % size of output field = # of particles

for m=1:length(x)  % assumed same as length(y)
        
    % get x,y as fractions of nx (remove periodic wrap-around)
    xl = mod(x(m)/dx, nx);
    yl = mod(y(m)/dy, ny);
        
    % get indeces of left/bottom grid point of cell containing
    % particle, min(i0,j0) = 1,1
    i0 = 1 + floor(xl);
    j0 = 1 + floor(yl);
        
    % get fractional position within cell, 0 <= ax,ay < 1
    ax = 1 + xl - i0;
    ay = 1 + yl - j0;
    
    wx = ones(2*(Iord+1),1); wy = wx;
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            if (i~=j) 
                wx(i+Iord+1) = wx(i+Iord+1)*(ax - j + bump)/(j - i);
                wy(i+Iord+1) = wy(i+Iord+1)*(ay - j + bump)/(j - i);
            end 
        end 
    end 
    
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            ig = 1 + mod( i0 + i - 1, nx );
            jg = 1 + mod( j0 + j - 1, nx );
            FI(m) = FI(m) + wx(i+Iord+1)*wy(j+Iord+1)*F(ig,jg);
        end 
    end     
end
      


