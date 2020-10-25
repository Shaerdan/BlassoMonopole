clc;
clear all;
global LogSwitch
LogSwitch = 1;
SourceNum = 2;
syms r theta psi theta0 psi0 I I1 I2 r1 r2 theta1 theta2 psi1 psi2
I = sym('I',[1,SourceNum]);
theta= sym('theta',[1,SourceNum]);
psi = sym('psi',[1,SourceNum]);
r = sym('r',[1,SourceNum]);
I_label = I;
cos_gamma = cos(theta).*cos(theta0)+sin(theta).*sin(theta0).*cos(psi-psi0);
dist = sqrt(1+r.^2-2*r.*cos_gamma);
phi_i  = (2./dist - log(1 - r.*cos_gamma + dist));
phi = sum(I.*phi_i);



% s = [r theta psi theta0];
Ntheta = 10; Npsi = 10;
QIncreSizeTheta = pi/Ntheta;
QIncreSizePsi = 2*pi/Npsi;
[theta0,psi0] =  meshgrid( 0:QIncreSizeTheta:pi,0:QIncreSizePsi:2*pi );
theta0 = theta0(:); psi0=psi0(:);
phi_i = subs(phi_i);
phi_norm = sqrt(sum(phi_i.^2,1));
phi_normalized = sum(subs(I.*phi_i)./phi_norm,2);
% norm_phi = ;
% phi_normalized = phi/norm_phi;
M = (Ntheta+1)^2;
dphi_dI = sym(zeros(M,SourceNum));
dphi_dr = sym(zeros(M,SourceNum));
dphi_dtheta = sym(zeros(M,SourceNum));
dphi_dpsi = sym(zeros(M,SourceNum));
for i = 1:SourceNum
    dphi_dI(:,i) = diff(phi_normalized,1,I(i));
    dphi_dr(:,i) = diff(phi_normalized,1,r(i));
    dphi_dtheta(:,i) = diff(phi_normalized,1,theta(i));
    dphi_dpsi(:,i)  = diff(phi_normalized,1,psi(i));
end

Jac = [dphi_dI dphi_dr dphi_dtheta dphi_dpsi];
%
I_exact = [1, -1];
r_exact = [0.5,0.5];
theta_exact = [0.9,0.901];
psi_exact = [0.5,0.501];

%
% I1 = 1; I2 = -1;
% r1 = 0.7; r2 = 0.7;
% theta1 = 0.5; theta2 = 0.501; psi1 = 0.5; psi2 = 0.501;
old = cell(4*SourceNum,1);
new = cell(4*SourceNum,1);
for i = 1:SourceNum
    new(i) = {I_exact(i)};
    new(i+SourceNum) =  {r_exact(i)};
    new(i+2*SourceNum) = {theta_exact(i)};
    new(i+3*SourceNum) = {psi_exact(i)};
    old(i) = {I(i)};
    old(i+SourceNum) =  {r(i)};
    old(i+2*SourceNum) = {theta(i)};
    old(i+3*SourceNum) = {psi(i)};
end

phi_normalized = subs(phi_normalized,old,new);

phi_normalized_num = double(phi_normalized);
surf(reshape(phi_normalized_num,[Ntheta+1 Npsi+1]))

Jac = double(subs(Jac,old,new));
pause(5)
Mesh.ThetaQ = theta0; Mesh.PsiQ=psi0; Lambda = 1e-2;SwitchNorm = 1;Measurement = phi_normalized_num;
X = [I_exact r_exact theta_exact psi_exact];
[f,JacNumeric] = GradPseudoLBFGSB(X,Measurement,Mesh,Lambda,SwitchNorm);
A = Jac - JacNumeric;
norm(A,inf)