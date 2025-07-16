function ETD_AC3D_characteristic_variable(N,M,K,space_interp)
format short e
format compact
N=128;
M=128;
L=128;
K=1000;
space_interp='linear';
%This program solve the Allen-Cahn equation 
%  u_t+ v*grad u = D(\Delta(u)-(u^3-u)/epsilon^2)
 
%Domain
%  [xb,xe] x [yb,ye] x [zb,ze]
%Number of grid sizes 
%   N x M x L
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
global z
global dt
global X
global Y
global Z
 
xb = 0; xe = 1; yb = 0; ye = 1; zb = 0; ze = 1;
T0 = 0; Te = 50; ord=2; bct=1;
 
hx = (xe-xb)/N;  hy = (ye-yb)/M; hz = (ze-zb)/M;
x = xb:hx:xe;  y = yb:hy:ye; z = zb:hz:ze; [X,Y,Z]=meshgrid(x,y,z);
dt = (Te-T0)/K;  T = T0:dt:Te;
fprintf(1,'\n *************************************************\n');
fprintf(1,'\n --- Parameters and Mesh Sizes ---\n');
fprintf(1,'\n epsilon = %e, D = %e\n',epsilon,D);
fprintf(1,'\n Nx = %d, Ny = %d, Nz = %d, Nt = %d\n',N,M,L,K);
fprintf(1,'\n hx = %e, hy = %e, hz = %e, dt = %e\n',hx,hy,hz,dt);
 
if bct==1
   fprintf(1,'\n Boundary Condition Type: Periodic\n');
end
if bct==2
   fprintf(1,'\n Boundary Condition Type: Neumann \n');
end
 
 
%Compute A,B,C
if bct==1
   NN = N; MM = M; LL = L;
   A = zeros(NN,NN);  B = zeros(MM,MM);  C=zeros(LL,LL); 
   for i=1:NN-1
       A(i,i) = -2;  A(i,i+1) = 1;  A(i+1,i) = 1;
   end
   A(NN,NN) = -2;  A(1,NN) = 1;  A(NN,1) = 1;
   for i=1:MM-1
       B(i,i) = -2;  B(i,i+1) = 1;  B(i+1,i) = 1;
   end
   B(MM,MM) = -2;  B(1,MM) = 1;  B(MM,1) = 1;
   for i=1:LL-1
       C(i,i) = -2;  C(i,i+1) = 1;  C(i+1,i) = 1;
   end
   C(LL,LL) = -2;  C(1,LL) = 1;  C(LL,1) = 1;
end
if bct==2
   NN = N+1; MM = M+1; LL = L+1;
   A = zeros(NN,NN);  B = zeros(MM,MM);  C=zeros(LL,LL);
   for i=1:NN-1
       A(i,i) = -2;  A(i,i+1) = 1;  A(i+1,i) = 1;
   end
   A(NN,NN) = -2;  A(1,2) = 2; A(NN,NN-1) = 2;
   for i=1:MM-1
       B(i,i) = -2;  B(i,i+1) = 1;  B(i+1,i) = 1;
   end
   B(MM,MM) = -2;  B(1,2) = 2;  B(MM,MM-1) = 2;
   for i=1:LL-1
       C(i,i) = -2;  C(i,i+1) = 1;  C(i+1,i) = 1;
   end
   C(LL,LL) = -2;  C(1,2) = 2;  C(LL,LL-1) = 2;
end

%Diagonalize A,B,C
A = D*A/(hx*hx);
[Px,Dx]=eig(A);
%Px; Dx;
Pxi = inv(Px);
B = D*B/(hy*hy);
[Py,Dy]=eig(B);
%Py; Dy;
Pyi = inv(Py);
C = D*C/(hz*hz);
[Pz,Dz]=eig(C);
%Pz; Dz;
Pzi = inv(Pz);
 
Het = zeros(NN,MM,LL); Het1 = zeros(NN,MM,LL);
phi1 = zeros(NN,MM,LL); phi2 = zeros(NN,MM,LL); 
for i=1:NN
    for j=1:MM
        for l=1:LL
            Hijl = Dx(i,i)+Dy(j,j)+Dz(l,l)-kappa;
            if abs(Hijl)>1.0e-3
               Het(i,j,l) = exp(Hijl*dt);
               Het1(i,j,l) = exp(Hijl*dt/2);
               phi1(i,j,l) = -(1-Het(i,j,l))/Hijl;
               phi2(i,j,l) = -(1-phi1(i,j,l)/dt)/Hijl;
            else
              Het(i,j,l) = 1;
              Het1(i,j,l)= 1;
              phi1(i,j,l) = dt;
              phi2(i,j,l) = dt/2;
           end
        end
    end
end
 
%Initialize U0
U = zeros(NN,MM,LL);
for i=1:NN
    for j=1:MM
        for l=1:LL
            U(i,j,l) = u0(x(i),y(j),z(l));
        end
    end
end

% U = 1.8*(rand(NN,MM,LL) - 0.5);
% U = load("U2.mat");
% U = struct2array(U);
F = zeros(NN,MM,LL,3);
if bct==1
V1 = zeros(NN+1,MM+1,LL+1,2);V2 = zeros(NN+1,MM+1,LL+1,2);V3=zeros(NN+1,MM+1,LL+1,2);
elseif bct==2
V1 = zeros(NN,MM,LL,2);V2 = zeros(NN,MM,LL,2);V3=zeros(NN,MM,LL,2);    
end
tic
 
%Evolve U along the time
 
TK = K;
vol =  0:1:TK;
egy =  0:1:TK;
maxi = 0:1:TK;
umin = 0:1:TK;
mass = 0:1:TK;
U_full = Recover(U,bct);
[vol(1),egy(1)] = comput_Vol_Egy(U_full,hx,hy,hz);
maxi(1)=max(abs(U_full(:)));
mass(1)=sum(U_full(:))*hx*hy*hz;
umin(1)=min(U_full(:));
for k=1:TK
    if ord==1 % ETD1-GODNOV splitting operater scheme (first-order scheme)
       [V1(:,:,:,1),V2(:,:,:,1),V3(:,:,:,1)] = Vec(X,Y,Z,k*dt);
       [Xc,Yc,Zc] = back_tracking(V1,V2,V3,1,bct);
       U_full1 = Recover(U,bct);
       U = interp3(x,y,z,U_full1,Xc,Yc,Zc,space_interp);
       
       F(:,:,:,1) = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,FF(U))));
       U_lin = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,U))).*Het;
       U = circle_z(Pz,circle_y(Py,circle_x(Px,U_lin-F(:,:,:,1)).*phi1));
    elseif ord==2 % ETD-RK2-Strang splitting operator scheme (second-order scheme)
    % B
    F1 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,FF(U))));
    U_lin1 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,U))).*Het;
    U = circle_z(Pz,circle_y(Py,circle_x(Px,U_lin1-F1.*phi1)));

    F2 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,FF(U))))-F1;
    U = U+circle_z(Pz,circle_y(Py,circle_x(Px,F2.*phi2)));
    % A
    [V1(:,:,:,1),V2(:,:,:,1),V3(:,:,:,1)] = Vec(X,Y,Z,k*dt);
    [Xc1,Yc1,Zc1] = back_tracking(V1,V2,V3,1,bct); % SSPRK1
    
    [Vx,Vy,Vz] = Vec(Xc1,Yc1,Zc1,k*dt-dt);
    V1(:,:,:,2) = Recover(Vx,bct); V2(:,:,:,2) = Recover(Vy,bct); V3(:,:,:,2) = Recover(Vz,bct);
    [Xc,Yc,Zc] = back_tracking(V1,V2,V3,2,bct);   % SSPRK2

    U_full = Recover(U,bct);
    U = interp3(x,y,z,U_full,Xc,Yc,Zc,space_interp);
    
    U = max(-beta,min(beta,U));

    % B
    F3 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,FF(U))));
    U_lin2 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,U))).*Het;
    U = circle_z(Pz,circle_y(Py,circle_x(Px,U_lin2-F3.*phi1)));

    F4 = circle_z(Pzi,circle_y(Pyi,circle_x(Pxi,FF(U))))-F3;
    U = U+circle_z(Pz,circle_y(Py,circle_x(Px, F4.*phi2)));
    end
    U_full = Recover(U,bct);
    mass(k+1) = sum(U_full(:))*hx*hy*hz;
    [vol(k+1),egy(k+1)] = comput_Vol_Egy(U_full,hx,hy,hz);
    maxi(k+1)=max(abs(U_full(:)));
    umin(k+1)=min(U_full(:));
end
 
wtime = toc;
fprintf (1,'\n MY_PROGRAM took %f seconds to run !!!\n',wtime);
 
fprintf(1,'\n Minimal value = %f, Maximal value = %f\n',min(min(min(U))),max(max(max(U))));
fprintf(1,'\n Initial energy = %e, Final energy is %e\n',egy(1),egy(TK+1));
fprintf(1,'\n Initial volume = %e, Final volume is %e\n',vol(1),vol(TK+1));
Radi = real((3/4*((mass-umin)./(maxi-umin))/pi).^(1/3));
 
%Draw the solution at the end time 
figure(1);
isosurface(x,y,z,U_full);
xlabel('X');
ylabel('Y');
zlabel('Z');
xticks(xb:0.2:xe);
yticks(yb:0.2:ye);
zticks(zb:0.1:ze);
axis([0 1 0 1 0 1]);
grid off;
set(figure(1), 'Position', [1000, 431, 438,389]); 
view([-51.3658919097676 24.4883244176014]);
% green_values = [0.8, 0.6, 0.4, 0.2, 0, 0.8, 0.6, 0.4, 0.2, 0.6, 0.4];
% green_colormap = [zeros(size(green_values))' green_values' zeros(size(green_values))'];
% colormap(green_colormap);
% % 构造 RGB 颜色映射矩阵（红=0，绿=自定义值，蓝=0）
% custom_green = [zeros(size(green_values,2))', green_values', zeros(size(green_values,2))'];
% % 应用颜色映射
% colormap(custom_green);
%view([0,0,1]);
 
%Draw the volume and energy curves along the time 
figure(2);
plot(T(1:TK+1),vol(1:TK+1),'.-');
xlabel('Time');
ylabel('Volume');

figure(3);
plot(T(1:TK+1),egy(1:TK+1),'.-');
xlabel('Time');
ylabel('Energy');

figure(4);
plot(T(1:TK+1),maxi(1:TK+1),'.-');
xlabel('Time');
ylabel('Maximum');

figure(5);
surf(x,y,U_full(:,:,round(LL/2)));
colormap jet;
shading interp;    
% axis off;
colorbar;

figure(6); 
surf(X(:,:,round(NN/2)), Y(:,:,round(MM/2)), Z(:,:,round(LL/2)), squeeze(U_full(:,:, round(LL/2))), 'EdgeColor', 'none'); 
colormap jet; 
grid off;set(figure(6), 'Position', [1000, 431, 438,389]); 
xlabel('X');
ylabel('Y');
zlabel('Z');
xticks(xb:0.2:xe);
yticks(yb:0.2:ye);
zticks(zb:0.1:ze);
axis([0 1 0 1 0 1]);
grid on;
view([-51.3658919097676 24.4883244176014]);
% ax = gca;ax.LineWidth = 1;

figure(7);
plot(T(1:TK+1),Radi(1:TK+1),'.-');
xlabel('Time');
ylabel('Radius');
end
 
function V = circle_x(P,U)
V = U;
for l=1:size(U,3)
    V(:,:,l) = P*U(:,:,l);
end
end

function V = circle_y(P,U)
V = U;
for i=1:size(U,1)
    U1(:,:) = U(i,:,:);
    V1 = P*U1;
    V(i,:,:) = V1;
end
end

function V = circle_z(P,U)
V = U;
for j=1:size(U,2)
    U1(:,:) = U(:,j,:);
    V1= U1*P'; 
    V(:,j,:) = V1;
end
end
  
function [vol, egy] = comput_Vol_Egy(U, hx, hy, hz)
global D poten; 
theta = 0.8; 
theta_c = 1.6;

[N, M, L] = size(U);
vol = 0;
egy1 = 0;  
egy2 = 0; 
egy3 = 0;

% 遍历三维网格单元
for i = 1:N-1
    for j = 1:M-1
        for k = 1:L-1
            % 计算当前立方体单元8个顶点的平均值
            uu = 0.125 * (U(i,j,k) + U(i+1,j,k) + U(i,j+1,k) + U(i+1,j+1,k) + ...
                        U(i,j,k+1) + U(i+1,j,k+1) + U(i,j+1,k+1) + U(i+1,j+1,k+1));
            
            % 体积积分项
            vol = vol + (uu + 1)/2;
            
            % 势能项
            egy1 = egy1 + (uu^2 - 1)^2;
            
            % 熵项
            egy3 = egy3 + theta/2 * ((1+uu)*log(uu+1) + (1-uu)*log(1-uu)) - theta_c/2 * uu^2;
            
            % 计算三个方向梯度（中心差分）
            udx = 0.25*( (U(i+1,j,k) - U(i,j,k)) / hx + ...
                        (U(i+1,j+1,k) - U(i,j+1,k)) / hx + ...
                        (U(i+1,j,k+1) - U(i,j,k+1)) / hx + ...
                        (U(i+1,j+1,k+1) - U(i,j+1,k+1)) / hx );
                    
            udy = 0.25*( (U(i,j+1,k) - U(i,j,k)) / hy + ...
                        (U(i+1,j+1,k) - U(i+1,j,k)) / hy + ...
                        (U(i,j+1,k+1) - U(i,j,k+1)) / hy + ...
                        (U(i+1,j+1,k+1) - U(i+1,j,k+1)) / hy );
                    
            udz = 0.25*( (U(i,j,k+1) - U(i,j,k)) / hz + ...
                        (U(i+1,j,k+1) - U(i+1,j,k)) / hz + ...
                        (U(i,j+1,k+1) - U(i,j+1,k)) / hz + ...
                        (U(i+1,j+1,k+1) - U(i+1,j+1,k)) / hz );
            
            % 梯度能量项
            egy2 = egy2 + (udx^2 + udy^2 + udz^2);
        end
    end
end

% 根据势能类型计算总能量
if poten == 1
    vol = vol * hx * hy * hz;
    egy1 = 0.25 * egy1;
    egy2 = D * 0.5 * egy2;
    egy = (egy1 + egy2) * hx * hy * hz;
elseif poten == 2
    vol = vol * hx * hy * hz;
    egy2 = D * 0.5 * egy2;
    egy = (egy3 + egy2) * hx * hy * hz;
end
end
 
function U_full = Recover(U,bct)
NN = size(U,1);  MM = size(U,2);  LL = size(U,3);
if bct==1
   U_full = zeros(NN+1,MM+1,LL+1); 
   U_full(1:NN,1:MM,1:LL) = U(1:NN,1:MM,1:LL);
   U_full(1:NN,1:MM,LL+1) = U_full(1:NN,1:MM,1);
   U_full(1:NN,MM+1,1:LL+1) = U_full(1:NN,1,1:LL+1);
   U_full(NN+1,1:MM+1,1:LL+1) = U_full(1,1:MM+1,1:LL+1);
end  
if bct==2
   U_full = U;
end
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
 
function uu = u0(x,y,z)
% uu = 0.9*sin(100*pi*x)*sin(100*pi*y)*sin(100*pi*z);
if (x-0.3)^2 + (y-0.3)^2 + (z-0.5)^2 <0.2^2
    uu = 0.9;
else
    uu = -0.9;
end
end
 
function [V1,V2,V3]=Vec(x,y,z,t)
V2=0.5-x;
V1=y-0.5;
V3=0*z;
% V1=exp(-t)*sin(2*pi*x);
% V2=exp(-t)*sin(2*pi*y);
% V3=exp(-t)*sin(2*pi*z);
end
 
function [Xc,Yc,Zc]=back_tracking(V1,V2,V3,ord,bct)
global dt
global X
global Y
global Z
if ord==1
 Xc=X-dt*V1(:,:,:,1);
 Yc=Y-dt*V2(:,:,:,1);
 Zc=Z-dt*V3(:,:,:,1);
elseif ord==2
 Xc=X-dt/2*(V1(:,:,:,1)+V1(:,:,:,2));
 Yc=Y-dt/2*(V2(:,:,:,1)+V2(:,:,:,2));
 Zc=Z-dt/2*(V3(:,:,:,1)+V3(:,:,:,2));
end
if bct==1
Xc = mod(Xc(1:end-1, 1:end-1, 1:end-1), 1); 
Yc = mod(Yc(1:end-1, 1:end-1, 1:end-1), 1);
Zc = mod(Zc(1:end-1, 1:end-1, 1:end-1), 1);
end
end