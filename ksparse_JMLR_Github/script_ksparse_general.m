%% ----------------------   SCRIPT CLUSTERING    --------------------------
fprintf('-----------------  SCRIPT CLUSTERING   ----------------------\n');
clear;
close all;
rng(1);

addpath(genpath('functions'));
addpath(genpath('data_prepared'));

%-------------Import data :
load('DataPatel.mat');
featuresnames = genenames;
X = X;
k = 5; 


%-------------Ksparse clustering :
param.LDA_ETA = 700;          % sparsity constraint
initialization = 'spectral';  % 'spectral' or 'PCA'
isTsne = 1;                   % 1 to display tsne, 0 otherwise
isSil = 1;                    % 1 to compute silhouette, 0 otherwise
% Advanced parameters :
param.h = k+8;
param.LOOP = 10;              % number of iteration
param.nb_kms = 40;            % number of kmeans-replicates
param.sigma = 150;            % used for initialisation with spectral
param.LDA_MAXITER = 50;       % be careful 10 with FISTA 
param.LDA_STEPSIZE = 1/1.001; % be careful 1/1.01 with Fista_AC

[Ysd,w,NormFrob] = ksparse(X,k,param,initialization,isTsne,isSil);
topGenes = top_genes_norm(w,featuresnames);

figure
plot(NormFrob,'LineWidth',2)
set(gca, 'YScale', 'log')
title('Frobenius norm')

figure
normGene = zeros(size(w,1),1);
for i = 1:size(w,1)
    normGene(i,1) = norm(w(i,:));
end
[normSort] = sort(normGene,'descend');
plot(normSort,'LineWidth',2)
title('sort( ||W(i,:)|| )')

    
%-------Eva Clusters 
Xw=X*w;
eva = evalclusters(Xw,'kmeans','Gap','KList',3:7);
display(eva);

%-------Size of each class (estimated) :
Ctab = {num2str(zeros(1,k))};
size_per_class = zeros(1,k);
for j = 1:k
    size_per_class(1,j) = sum(Ysd==j);
    tmp = ['C' int2str(j)];
    Ctab{j} = tmp;
end
size_per_class = mat2dataset(size_per_class,'VarNames',Ctab,'ObsNames',{'ksparse'});
display(size_per_class);


