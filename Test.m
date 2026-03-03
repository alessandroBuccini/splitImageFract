%% Clear workspace
clc
clear
close all

%% Choose example
choice=menu('Please choose one of the following options:','Circumference', 'Triangle', 'Writings', 'DomusDeJanas');
test_case_list = {'Circumference', 'Triangle', 'Writings', 'DomusDeJanas'};
test_case = test_case_list{choice};

%% Construct example
load([test_case,'.mat'])
if exist('G_true','var')
    D=createData(S_true,G_true,sigma,nu);
    B_true=D-G_true;
end

%% Set parameters
options.waitbar='on';
options.tol=1e-6;
options.iter=200;

%% Recover
l_mu=length(mu);
BB=cell(l_mu,1);
GG=cell(l_mu,1);
OUTINFO=cell(l_mu,1);
ERR=zeros(l_mu,2);

for j=1:length(mu)
    [B,G,outInfo]=splitImageFract(D,mu(j),options);
    BB{j}=B;
    GG{j}=G;
    OUTINFO{j}=outInfo;
    if exist('G_true','var')
    for p=1:2
        ERR(j,p)=norm(G(:)/norm(G(:),p)-G_true(:)/norm(G_true(:),p),p);
    end
    end
end

%% Show results
if ~exist('CLIM','var')
    CLIM=[];
end
figure(1),
imshow(D,[]), title('Archaelogical Surface'), colormap hsv(256);
figure(2),
if choice==4 
    tr=25;
        for j=1:length(mu)
        subplot(2,length(mu),j), imshow(BB{j}(tr:end-tr,tr:end-tr),CLIM), title(['Recovered B, \mu = ', num2str(mu(j),2)]), xlabel(['iterations = ',num2str(OUTINFO{j}.iter)])
        subplot(2,length(mu),length(mu)+j), imshow(GG{j}(tr:end-tr,tr:end-30),CLIM), title(['Recovered G, \mu = ', num2str(mu(j),2)]),xlabel(['||D-G-B|| = ',num2str(OUTINFO{j}.cstr(end),2)]),        
        end

else
    subplot(2,length(mu)+1,1), imshow(D-G_true,CLIM), title('Exact B'), colormap hsv(256);
    subplot(2,length(mu)+1,length(mu)+2), imshow(G_true,CLIM), title('Exact G'), colormap hsv(256);
    for j=2:length(mu)+1
        subplot(2,length(mu)+1,j), imshow(BB{j-1},CLIM), title(['Recovered B, \mu = ', num2str(mu(j-1),2)]), xlabel(['iterations = ',num2str(OUTINFO{j-1}.iter)])
        subplot(2,length(mu)+1,length(mu)+1+j), imshow(GG{j-1},CLIM), title(['Recovered G, \mu = ', num2str(mu(j-1),2)]),xlabel(['||D-G-B|| = ',num2str(OUTINFO{j-1}.cstr(end),2)])
    end
end
if choice==4
    colormap(flipud(hsv(256)));
else
    colormap hsv(256)
end

%% Show Errors
if exist('G_true','var')
    [~,i_best]=min(ERR);
    figure(3)
    loglog(mu,ERR,mu(i_best(1)),ERR(i_best(1),1),'b*',mu(i_best(2)),ERR(i_best(2),2),'r*')
    xlabel('$\mu$',Fontsize=15,Interpreter='latex')
    xlim([mu(1),mu(end)])
    legend('$f_1(G)$','$f_2(G)$','Location','best',Fontsize=15,Interpreter='latex')
end