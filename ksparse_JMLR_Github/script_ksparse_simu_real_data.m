%% ----------------------   SCRIPT CLUSTERING    --------------------------
fprintf('-----------------  SCRIPT CLUSTERING   ----------------------\n');
clear;
close all;
rng(1);

addpath(genpath('functions'));
addpath(genpath('Data_prepared'));

%-------------Import real data :

%load('DataPatel.mat');
% load('DataUsoskin.mat');
%load('DataKlein.mat');
%load('DataZeisel9.mat');
% 
% k = nb_clusters;
% featuresnames = genenames;
% m = size(X,1);
% d = size(X,2);


%-------------Import simu data :

load('data_simu_10000_k4.mat');
%load('data_simu_15000_k4.mat');
k = ksimu;
X = Xcell;
featuresnames = genenames;
m = size(X,1);
d = size(X,2);



%___________________________CLUSTERING : 
%-------------PCA_kmeans :
fprintf('-> PCA_kmeans... \n');
nb_kms = 10;
nb_components_pca = 5;
tic
[proj_axes,Ykms] = kmeans_PCA(X,k,nb_components_pca,nb_kms);
timePCA = toc;
Ykms = match_names(YR,Ykms,k);

%-------------Spectral_sym : 
fprintf('-> Spectral_sym... \n');
sigma = 150;
tic
[Ysym,projSpc] = spcl_normalized_sym(X,k,sigma,nb_kms);
timeSpectral = toc;
Ysym = match_names(YR,Ysym,k);

%------------- SIMLR :
fprintf('-> SIMLR... \n');
tic
[Ysimlr, S, F, XwSIMLR, alpha] = SIMLR(X,k,10);
timeSIMLR = toc;
Ysimlr = match_names(YR,Ysimlr,k);


%-------------Ksparse :
% Patel   : k=5, ETA = 5000,   
% Usoskin : k=4, ETA = 5000,  
% Klein   : k=4, ETA = 25000, 
% Zeisel  : k=9, ETA = 13000, 
param.LDA_ETA = 5000;          % sparsity constraint
isTsne = 0;                    % 1 to display tsne, 0 otherwise
param.LOOP = 10;               % number of loops
param.h = k+4;                 % Size of the centroids
param.nb_kms = 40;             % number of kmeans-replicates
param.sigma = 150;             % used for initialisation with spectral
initialization = 'spectral';   % 'spectral' or 'PCA'
param.LDA_MAXITER = 50;        % be careful 10 with FISTA 
param.LDA_STEPSIZE = 1/1.001;  % be careful 1/1.01 with Fista_AC  

tic
[Ysd,w,NormFrob,Si_sd,XW] = ksparse(X,k,param,initialization,isTsne);
timeSparse = toc;
Ysd = match_names(YR,Ysd,k);
topGenes = top_genes_norm(w,featuresnames);

figure
plot(NormFrob,'LineWidth',2)
set(gca, 'YScale', 'log')
title('Frobenius norm')

 
        
    
%_____________________________RESULTS  :
%-------Size of each class (estimated) :
Ctab = {num2str(zeros(1,k))};
size_class = zeros(5,k);
for j = 1:k
    size_class(1,j) = sum(YR==j);
    size_class(2,j) = sum(Ykms==j);
    size_class(3,j) = sum(Ysym==j);
    size_class(4,j) = sum(Ysimlr==j);
    size_class(5,j) = sum(Ysd==j);
    tmp = ['C' int2str(j)];
    Ctab{j} = tmp;
end
size_class = mat2dataset(size_class,'VarNames',Ctab,'ObsNames',{'Reality','PCA_kmeans','Spectral_sym','SIMLR','Ksparse'});
display(size_class);

%-------ARI and Accuracy :
Tacc = zeros(4,3);
Tacc(1,1) = RandIndex(YR,Ykms);
Tacc(1,2) = nmi(YR,Ykms);
Tacc(1,3) = compute_accuracy(YR,Ykms,k);
Tacc(2,1) = RandIndex(YR,Ysym);
Tacc(2,2) = nmi(YR,Ysym);
Tacc(2,3) = compute_accuracy(YR,Ysym,k);
Tacc(3,1) = RandIndex(YR,Ysimlr);
Tacc(3,2) = nmi(YR,Ysimlr);
Tacc(3,3) = compute_accuracy(YR,Ysimlr,k);
Tacc(4,1) = RandIndex(YR,Ysd);
Tacc(4,2) = nmi(YR,Ysd);
Tacc(4,3) = compute_accuracy(YR,Ysd,k);
acctab = {num2str(zeros(1,k+2))}; acctab{1} = 'ARI'; acctab{2} = 'NMI'; acctab{3} = 'global_accuracy';
disp_accuracy = mat2dataset(Tacc,'VarNames',acctab,'ObsNames',{'PCA_kmeans','Spectral_sym','SIMLR','Ksparse'});
display(disp_accuracy);

%-------Time :
tabTime = [timePCA; timeSpectral; timeSIMLR; timeSparse];
tabTime = mat2dataset(tabTime,'VarNames',{'time'},'ObsNames',{'PCA_kmeans','Spectral_sym','SIMLR','Ksparse'});
display(tabTime);

%-------Silhouette :
fprintf('-> Computing silhouette... \n');
figure('name','Silhouettes')
subplot(1,4,1)
[Si_PCA,h1] = silhouette(proj_axes,Ykms);
title('PCA kmeans')
subplot(1,4,2)
[Si_Spec,h2] = silhouette(projSpc,Ysym);
title('Spectral sym')
subplot(1,4,3)
[Si_Simlr,h3] = silhouette(XwSIMLR,Ysimlr);
title('SIMLR')
subplot(1,4,4)
[Si_sd,h4] = silhouette(XW,Ysd);
title('Ksparse')
Mean_sil = [mean(Si_PCA); mean(Si_Spec); mean(Si_Simlr); mean(Si_sd)];
Mean_sil = mat2dataset(Mean_sil,'VarNames',{'Silhouette_average'},'ObsNames',{'PCA_kmeans','Spectral_sym','SIMLR','ksparse'});
display(Mean_sil);


%-------TSNE graph :
Ykms2 = Ykms;
Ysym2 = Ysym;
Ysd2 = Ysd;
Ysimlr2 = Ysimlr;
Ykms2(Ykms2~=YR) = 0;
Ysym2(Ysym2~=YR) = 0;
Ysimlr2(Ysimlr2~=YR) = 0;
Ysd2(Ysd2~=YR) = 0;
Ykms2 = Ykms2 + ones(size(YR));
Ysym2 = Ysym2 + ones(size(YR));
Ysimlr2 = Ysimlr2 + ones(size(YR));
Ysd2 = Ysd2 + ones(size(YR));
fprintf('-> Computing tsne figures... \n');
figure('name','TSNE graph')
subplot(1,4,1)
T = tsne(proj_axes);
gscatter(T(:,1),T(:,2),Ykms2)
title_PCA = sprintf('PCA kmeans %d errors (red)', sum(Ykms2==1));
title(title_PCA);
subplot(1,4,2)
T = tsne(projSpc);
gscatter(T(:,1),T(:,2),Ysym2)
title_spectral = sprintf('Spectral sym %d errors (red)', sum(Ysym2==1));
title(title_spectral);
subplot(1,4,3)
T = tsne(XwSIMLR);
gscatter(T(:,1),T(:,2),Ysym2)
title_simlr = sprintf('SIMLR %d errors (red)', sum(Ysym2==1));
title(title_simlr);
subplot(1,4,4)
T = tsne(XW);
gscatter(T(:,1),T(:,2),Ysd2)
title_sparse = sprintf('Ksparse %d errors (red)', sum(Ysd2==1));
title(title_sparse);

