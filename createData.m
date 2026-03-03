function D = createData(S,G,sigma,nu)
% This functions create a synthetic dataset for glyph extraction. The
% archaelogical surface is obtained as the sum of a smooth part, denoted by
% xs, and a sparse part, denoted by xg. Roughness component is added to 
% simulate more realistic data. The smoothness level of roughness is denoted
% by sigma and its level by nu.
%
% x=createData(xs,xg,sigma,nu);
%
% INPUT:
% S: smooth part of the archaelogical surface (must be of the same size as
%     G);
% G: sparse part of the archaelogical surface (must be of the same size as
%     S);
% sigma: roughness smoothness;
% nu: roughness level. The roughness affects only the area where abs(G)>eps.
%                      and near pixels; 
%
% OUTPUT:
% D: archaelogical surface.
%
% This function is associated to the paper:
% A. Azzarelli and A.Buccini, ...


% Check input
[nS,mS]=size(S);
[nG,mG]=size(G);
if nS~=nG || mS~=mG
    error('The dimensions of B and G should match');
end

% Create roughness component
mask=fspecial('gaussian',fix(5*sigma),sigma);
R=randn(nS,mS);
err=conv2(R,mask,'same');
err=err/norm(err,'fro');
if min(nG,mG)<20
    glyphzone = abs(G)>eps;
else
    glyphzone = 1-conv2(1-(abs(G)<=eps),fspecial('gaussian',10,2),'same');
end
% Combine the elements
B=S+nu*norm(S+G,'fro')*err.*glyphzone;
D=B+G;