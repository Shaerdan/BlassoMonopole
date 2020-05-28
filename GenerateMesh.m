function [Mesh] = GenerateMesh(MeshSize_R,MeshSize_Theta,MeshSize_Psi)
%UNTITLED Summary of this function goes here
%   Generate Mesh
%   Mesh.ThetaQ:  a 2d matrix meshing the theta commponent of the point Q \in R^2;
%   Mesh.PsiQ:    a 2d matrix meshing the psi commponent of the point Q \in R^2;
%   Mesh.RadiusS:   a 3d array(tensor) meshing of the radius commponent of the point S \in R^3;
%   Mesh.ThetaS:    a 3d array(tensor) meshing of the theta commponent of the point S \in R^3;
%   Mesh.PsiS:    a 3d array(tensor) meshing of psi commponent of the point S \in R^3;

hr = 1/MeshSize_R;
ht = pi/MeshSize_Theta;
hp = 2*pi/MeshSize_Psi;

r   = 0:hr:0.999; 
Mesh.ThetaLine=0:ht:pi; 
Mesh.PsiLine = 0:hp:2*pi; 

[Mesh.ThetaQ,Mesh.PsiQ] = meshgrid(Mesh.ThetaLine,Mesh.PsiLine);
Mesh.RadiusS = r; Mesh.ThetaS = Mesh.ThetaLine; Mesh.PsiS = Mesh.PsiLine;
Mesh.LengthRadiusS = length(Mesh.RadiusS); 
Mesh.LengthThetaS = length(Mesh.ThetaS); 
Mesh.LengthPsiS = length(Mesh.PsiS);

end

