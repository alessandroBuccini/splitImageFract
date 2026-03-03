function [B,G,outInfo]=splitImageFract(D,mu,options)
% This function splits the surface D into a smooth component B and a
% sparse component G solving the minimization problem
%
% (B,G)= argmin ||L^alpha*b||_2^2+mu*||g||_1 subj. to D=B+G,
%
% where L is the Laplacian operator, b=vec(B) and g=vec(G)
% using the ADMM algorithm.
%
% [B,G]=splitImageFract(D,mu,options);
%
% INPUT:
%      D: surface to split;
%      mu: sparsity parameter;
%      options: structure containing the parameters for running the
%               algorithm (optional)
%             .alpha: fractional exponent of the Laplacian (if empty or not
%                     populated an heuristic rule is used for its estimation);
%             .rho: augmentation parameter for the ADMM (default rho=1);
%             .iter: maximum number of iterations (default iter=100);
%             .tol: tollerance of the stopping criterion, set tol=0 to
%                   force the method to perform all the iterations (default
%                   tol=1e-4);
%             .waitbar: flag to enable a waitbar. If options.waitbar=='on',
%                       a waitbar is displayed (default waitbar='off').
%
% OUTPUT
%      B: smooth component;
%      G: sparse component;
%      outInfo: structure containing informations on the iterations
%             .alpha: value of alpha
%             .J: value of the minimized functional at each iteration;
%             .cstr: value of the residual at each iteration (D-B-G);
%             .iter: number of performed iterations.
%
% This function is associated to the paper:
% A. Azzarelli and A.Buccini, ...

% Check input arguments
if nargin<2
    error('Not enough input arguments')
end
if nargin==2
    options.alpha=determineAlpha(D);
    options.rho=1;
    options.iter=100;
    options.tol=1e-4;
    options.waitbar='off';
end

if ~isscalar(mu)
    mu=1;
    warning(['The value of mu was not acceptable. It has been set to ',num2str(mu)]);
end

% Extract informations from the options structure
% alpha
if ~isfield(options,'alpha') || isempty(options.alpha)
    [alpha,L]=determineAlpha(D);
else
    alpha=options.alpha;
end
if ~isnumeric(alpha)
    [alpha,L]=determineAlpha(D);
    warning(['The value of options.alpha was not acceptable. It has been set to ',num2str(alpha)]);
end
% rho
if ~isfield(options,'rho') || isempty(options.rho)
    rho=1;
else
    rho=options.rho;
end
if ~isnumeric(rho) || rho<=0
    rho=1;
    warning('The value of options.rho was not acceptable. It has been set to 1');
end
% iter
if ~isfield(options,'iter') || isempty(options.iter)
    iter=100;
else
    iter=options.iter;
end
if ~isnumeric(iter) || iter<1
    iter=100;
    warning('The value of options.iter was not acceptable. It has been set to 100.');
end
% tol
if ~isfield(options,'tol') || isempty(options.tol)
    tol=1e-4;
else
    tol=options.tol;
end
if ~isnumeric(tol) || tol<0 || tol>1
    tol=1e-4;
    warning('The value of options.tol was not acceptable. It has been set to 1e-4');
end
% waitbar
try
    flagW=double(strcmp(options.waitbar,'on'));
catch
    flagW=0;
end

% Create output structure
if nargout==3
    outInfo.alpha=alpha;
    outInfo.J=zeros(iter,1);
    outInfo.cstr=zeros(iter,1);
    outInfo.iter=iter;
    flagOut=1;
else
    flagOut=0;
end

% Set initial guesses for ADMM
[n,m]=size(D);
B=D;
G=zeros(n,m);
Y=zeros(n,m);

% Compute the eigenvalues of L^alpha
if ~exist('L','var')
    L=[0 -1 0; -1 4 -1; 0 -1 0]; % stencil for the two-dimensional Laplacian operator
    L=padarray(L,[n,m]-size(L),0,'post'); % padding to make the stencil as big as the surface
    L=fft2(circshift(L,1-[2,2]));% Computation of the eigenvalues of the Laplacian matrix
end
L=L.^alpha;

% Show waitbar
if flagW
    h=waitbar(0,['splitImageFract - mu = ',num2str(mu,2),' - Computations in progress...']);
end

% Perform the ADMM iterations
for k=1:iter
    % Store the previous iterations for the stopping criterion
    B_old=B;
    G_old=G;

    % Solve the two subproblems
    [B,B_hat]=solve_B(D,G,Y,rho,L);
    G=solve_G(D,B,mu,Y,rho);

    % Update of the Lagrangian multiplier
    Y=Y+rho*(B+G-D);

    % Fill the output structure
    if flagOut
        outInfo.J(k)=0.5*norm(L.*B_hat,'fro')^2/numel(B_hat)+mu*norm(G(:),1);
        outInfo.cstr(k)=norm(B+G-D,'fro');
    end

    % Update waitbar
    if flagW
        waitbar(k/(iter-1));
    end

    % Check the stopping criterion
    if norm(B_old-B,'fro')<tol*norm(B_old,'fro') && norm(G_old-G,'fro')<tol*norm(G_old,'fro')
        if flagOut
            outInfo.J=outInfo.J(1:k);
            outInfo.cstr=outInfo.cstr(1:k);
            outInfo.iter=k;
        end
        break
    end
end

% Ensure that the results are real
B=real(B);
G=real(G);

% Close waitbar
if flagW
    try
        close(h)
    catch
        warning('Waitbar not closed.')
    end
end
end

function [B,B_hat]=solve_B(D,G,Y,rho,L)
rhs=fft2(rho*(D-G)-Y);
B_hat=rhs./(abs(L).^2+rho);
B=ifft2(B_hat);
end

function G=solve_G(D,B,mu,Y,rho)
Z=D-B-Y/rho;
G=abs(Z)-mu/rho;
G=max(G,0);
G=sign(Z).*G;
end

function [alpha,L]=determineAlpha(D)
[n,m]=size(D);
L=[0 -1 0; -1 4 -1; 0 -1 0];
L=padarray(L,[n,m]-size(L),0,'post');
L=fft2(circshift(L,1-[2,2]));
xhat=fft2(D);
alpha=fminbnd(@(alpha) norm(L.^alpha.*xhat)^2,0,2);
end