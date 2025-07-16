function U_full=ETD_AC2D_characteristic_variable(N,M,K,space_interp)
format short e
format compact
N=512;
M=512;
K=80;
space_interp='spline';
%This program solve the Allen-Cahn equation 
%  u_t+ v*grad u = D(\Delta(u)-(u^3-u)/epsilon^2)
 
%Domain
%  [xb,xe] x [yb,ye]
%Number of grid sizes 
%   N x M 
%Time interval
%  [T0,Te] 
%Number of step sizes
%   K
 
%Allen-Cahn parameters
global epsilon;
global D;
global kappa;
global theta;
global theta_c;
global poten;
global beta;

epsilon = 0.01;
D = epsilon^2;
theta = 0.8; 
theta_c = 2*theta; 
poten=2;

if poten==1
    kappa = 2;
    beta = 1;
elseif poten==2
    kappa = 8.02;
    beta = 0.9575;
end
 
global x
global y
global dt
global X
global Y
 
xb = 0; xe = 1; yb = 0; ye = 1; T0 = 0; Te = 8; ord=2;
 
hx = (xe-xb)/N;  hy = (ye-yb)/M;
x = xb:hx:xe;  y = yb:hy:ye; [X,Y]=meshgrid(x,y);
dt = (Te-T0)/K;  T = T0:dt:Te;
fprintf(1,'\n *************************************************\n');
fprintf(1,'\n --- Parameters and Mesh Sizes ---\n');
fprintf(1,'\n epsilon = %e, D = %e\n',epsilon,D);
fprintf(1,'\n Nx = %d, Ny = %d, Nt = %d\n',N,M,K);
fprintf(1,'\n hx = %e, hy = %e, dt = %e\n',hx,hy,dt);
 
fprintf(1,'\n Boundary Condition Type: Periodic\n');
 
 
%Compute A and B
 
   NN = N; MM = M;
   A = zeros(NN,NN);  B = zeros(MM,MM);
   for i=1:NN-1
       A(i,i) = -2;  A(i,i+1) = 1;  A(i+1,i) = 1;
   end
   A(NN,NN) = -2;  A(1,NN) = 1;  A(NN,1) = 1;
   for i=1:MM-1
       B(i,i) = -2;  B(i,i+1) = 1;  B(i+1,i) = 1;
   end
   B(MM,MM) = -2;  B(1,MM) = 1;  B(MM,1) = 1;
 
%Diagonalize A,B
A = D*A/(hx*hx);
[Px,Dx]=eig(A);
%Px; Dx;
Pxi = inv(Px);
B = D*B/(hy*hy);
[Py,Dy]=eig(B);
%Py; Dy;
Pyi = inv(Py);
 
Het = zeros(NN,MM); Het1 = zeros(NN,MM);
phi = zeros(NN,MM,3); 
for i=1:NN
    for j=1:MM      
        Hij = Dx(i,i)+Dy(j,j)-kappa;
        if abs(Hij)>1.0e-3
            Het(i,j) = exp(Hij*dt);
            Het1(i,j) = exp(Hij*dt/2);
            phi(i,j,1) = -(1-Het(i,j))/Hij;
            phi(i,j,2) = -(1-Het1(i,j))/Hij;
            phi(i,j,3) = -(1-phi(i,j,1)/dt)/Hij;
        else
            Het(i,j) = 1;
            Het1(i,j)= 1;
            phi(i,j,1) = dt;
            phi(i,j,2) = dt/2;
            phi(i,j,3) = dt/2;
        end
    end
end
 
%Initialize U0
U = zeros(NN,MM);
for i=1:NN
    for j=1:MM
        U(i,j) = u0(x(i),y(j));
    end
end

% U = 1.8*rand(NN,MM) - 0.9;
% U = load("U1.mat");
% U = struct2array(U);
F = zeros(NN,MM,3);
V1 = zeros(NN+1,MM+1,2);V2 = zeros(NN+1,MM+1,2);
tic
 
%Evolve U along the time
 
TK = K;
vol =  0:1:TK;
egy =  0:1:TK;
maxi = 0:1:TK;
umin = 0:1:TK;
mass = 0:1:TK;
U_full = Recover(U);
[vol(1),egy(1)] = comput_Vol_Egy(U_full,hx,hy);
maxi(1)=max(abs(U_full(:)));
umin(1)=min(U_full(:));
mass(1)=sum(U_full(:))*hx*hy;
for k=1:TK
    if ord==1 % ETD1-GODNOV splitting operater scheme (first-order scheme)
       [V1(:,:,1),V2(:,:,1)] = Vec(X,Y,k*dt);
       [Xc,Yc] = back_tracking(V1,V2,1);
       U_full1 = Recover(U);
       U = interp2(x,y,U_full1,Xc,Yc,space_interp);
       
       F(:,:,1) = circle_y(Pyi,circle_x(Pxi,FF(U)));
       U_lin = circle_y(Pyi,circle_x(Pxi,U)).*Het;
       U = circle_y(Py,circle_x(Px,U_lin-F(:,:,1).*phi(:,:,1)));
    elseif ord==2 % ETD-RK2-Strang splitting operator scheme (second-order scheme)
    % B
    F1 = circle_y(Pyi,circle_x(Pxi,FF(U)));
    U_lin1 = circle_y(Pyi,circle_x(Pxi,U)).*Het;
    U = circle_y(Py,circle_x(Px,U_lin1-F1.*phi(:,:,1)));

    F2 = circle_y(Pyi,circle_x(Pxi,FF(U)))-F1;
    U = U+circle_y(Py,circle_x(Px,F2.*phi(:,:,3)));
    % A
    [V1(:,:,1),V2(:,:,1)] = Vec(X,Y,k*dt);
    [Xc1,Yc1] = back_tracking(V1,V2,1); % SSPRK1

    % U_full = Recover(U);
    % U = interp2(x,y,U_full,Xc1,Yc1,space_interp);

    [Vx,Vy] = Vec(Xc1,Yc1,k*dt-dt); 
    V1(:,:,2) = Recover(Vx); V2(:,:,2) = Recover(Vy);
    [Xc,Yc] = back_tracking(V1,V2,2);   % SSPRK2

   
    U_full = Recover(U);
    U = interp2(x,y,U_full,Xc,Yc,space_interp);

    
    U = max(-beta,min(beta,U));

    % B
    F3 = circle_y(Pyi,circle_x(Pxi,FF(U)));
    U_lin2 = circle_y(Pyi,circle_x(Pxi,U)).*Het;
    U = circle_y(Py,circle_x(Px,U_lin2-F3.*phi(:,:,1)));

    F4 = circle_y(Pyi,circle_x(Pxi,FF(U)))-F3;
    U = U+circle_y(Py,circle_x(Px, F4.*phi(:,:,3)));
    end
    U_full = Recover(U);
    [vol(k+1),egy(k+1)] = comput_Vol_Egy(U_full,hx,hy);
    maxi(k+1)=max(abs(U_full(:)));
    umin(k+1)=min(U_full(:));
    mass(k+1)=sum(U_full(:))*hx*hy;
end
 
wtime = toc;
fprintf (1,'\n MY_PROGRAM took %f seconds to run !!!\n',wtime);
 
fprintf(1,'\n Minimal value = %f, Maximal value = %f\n',min(min(U)),max(max(U)));
fprintf(1,'\n Initial energy = %e, Final energy is %e\n',egy(1),egy(TK+1));
fprintf(1,'\n Initial volume = %e, Final volume is %e\n',vol(1),vol(TK+1));
Radi = real(((mass-umin)./(maxi-umin))/pi).^(1/2);
 
%Draw the solution at the end time 
fig1 = figure(1);
surf(x,y,U_full','linestyle','none');
xlabel('X');
ylabel('Y');
axis equal;
colorbar;
%view([0,0,1]);
 
%Draw the volume and energy curves along the time 
fig2 = figure(2);
plot(T(1:TK+1),vol(1:TK+1),'.-');
xlabel('Time');
ylabel('Volume');
fig3 = figure(3);
plot(T(1:TK+1),egy(1:TK+1),'.-');
xlabel('Time');
ylabel('Energy');
fig4 = figure(4);
plot(T(1:TK+1),maxi(1:TK+1),'.-');
xlabel('Time');
ylabel('Maximum');

fig5 = figure(5);
pcolor(x,y,U_full');
colormap jet;
shading interp;    
axis off;
colorbar;

figure(7);
plot(T(1:TK+1),Radi(1:TK+1),'.-');
xlabel('Time');
ylabel('Radius');
end
 
function V = circle_x(P,U)
V = P*U;
end
 
function V = circle_y(P,U)
V = U*P';
end
  
function [vol,egy] = comput_Vol_Egy(U,hx,hy)
global D poten; theta = 0.8; theta_c = 1.6;
N = size(U,1);  M = size(U,2);
vol = 0;
egy1 = 0;  egy2 = 0; egy3 = 0;
for i=1:N-1
    for j=1:M-1
        uu = 0.25*(U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1));
        vol = vol+(uu+1)/2;
        egy1 = egy1+(uu*uu-1)^2;
        egy3 = egy3 + theta/2*((1+uu)*log(uu+1) + (1-uu)*log(1-uu)) - theta_c/2*(uu)^2;
        udx = 0.5*(U(i+1,j)-U(i,j)+U(i+1,j+1)-U(i,j+1))/hx;
        udy = 0.5*(U(i,j+1)-U(i,j)+U(i+1,j+1)-U(i+1,j))/hy;
        egy2 = egy2+(udx*udx+udy*udy);
    end
end
if poten==1
vol = vol*hx*hy;
egy1 = 0.25*egy1;
egy2 = D*0.5*egy2;
egy = (egy1+egy2)*hx*hy;
end
if poten==2
vol = vol*hx*hy;
egy2 = D*0.5*egy2;
egy = (egy3+egy2)*hx*hy;
end
end
 
function U_full = Recover(U)
NN = size(U,1);  MM = size(U,2);
 
U_full = zeros(NN+1,MM+1);
U_full(1:NN,1:MM) = U(1:NN,1:MM);
U_full(1:NN,MM+1) = U_full(1:NN,1);
U_full(NN+1,1:MM+1) = U_full(1,1:MM+1);
end
 
function V = FF_org(U)
global theta_c poten theta;
if poten==1
V = U.*(U.*U-1); 
end
if poten==2
V = -theta/2*(log(1-U)-log(1+U))-theta_c*U;   
end
end

function V = FF(U)
global kappa;
V = FF_org(U) - kappa*U;
end
 
function uu = u0(x,y)
% uu = 0.9*sin(100*pi*x)*sin(100*pi*y);
if (x-0.3)^2 + (y-0.3)^2  <0.2^2
    uu = 0.9;
else
    uu = -0.9;
end

end
 
function [V1,V2]=Vec(x,y,t)
V2=0.5-x;
V1=y-0.5;
% V1=exp(-t)*sin(2*pi*x);
% V2=exp(-t)*sin(2*pi*y);
end
 
function [Xc,Yc]=back_tracking(V1,V2,ord)
global dt
global X
global Y
if ord==1
 Xc=X-dt*V1(:,:,1);
 Yc=Y-dt*V2(:,:,1);
elseif ord==2
 Xc=X-dt/2*(V1(:,:,1)+V1(:,:,2));
 Yc=Y-dt/2*(V2(:,:,1)+V2(:,:,2));
end
Xc = mod(Xc(1:end-1,1:end-1),1);
Yc = mod(Yc(1:end-1,1:end-1),1);
end