function [Y,w,NormFrob,Si_sd,XW] = ksparse(X,k,param,initialization,isTsne)
%----------------- INPUT
% X                  : data matrix
% k                  : number of class
% param.LDA_ETA      : sparsity constraint
% param.LOOP         : number of iteration
% param.LDA_STEPSIZE : gradient stepsize
% param.LDA_MAXITER  : number of iteration in LDA_Multiclass
% param.nb_kms       : number of kmeans-replicates
% param.sigma        : used for initialisation with spectral
% initialisation     : 'spectral' or 'PCA'
% isTsne             : 1 to plot tsne, 0 otherwise
% isSil              : 1 to plot silhouette, 0 otherwise
%----------------- OUTPUT 
% Y                  : labels
% w                  : Sparse LDA directions for each class
% NormFrob           : Frobenius norm evolution of (Ystar-XW)
%__________________________________________________________________________

fprintf('-> ksparse clustering...\n');
                          
n = size(X,1);
d = size(X,2);
X = X/norm(X);            % normalization => beta = 1
Pi = eye(n)-1/n*ones(n,1)*ones(n,1)';
X = Pi*X;                 % Projection Vector

initSize = size(X,2);     % number of gene before sparsity
NormFrob = zeros(1,param.LOOP+1);


%------- Initialization :
if strcmp(initialization,'PCA')
    [~,score] = pca(X,'NumComponents',param.h);
    [Yk,mu] = kmeans(score,k,'Replicates',param.nb_kms);
end
if strcmp(initialization,'spectral')
    [Yk,~,mu] = spcl_normalized_sym_h(X,k,param.sigma,param.nb_kms,param.h);
end

%mu = eye(k,param.h);

%------- Sparse Discriminative :
Winit = ones(d,param.h);
%Sinit=param.LDA_ETA/(d*param.h);
Sinit = 1.0/(d*param.h);
Winit=Winit*Sinit;

XtX = X'*X;
Y = class2indicator(Yk,k);
y_star = Y*mu;


NormFrob(1) = trace((y_star-X*Winit)*(y_star-X*Winit)');
fprintf('processing loop : ');
for i = 1:param.LOOP
    fprintf(' %d ',i);
    Xty = X'*y_star;
    %w = LDA_FISTA_AC_MEX(k,param.LDA_ETA,param.LDA_STEPSIZE,param.LDA_MAXITER,XtX,Xty,Winit);
    w = LDA_FISTA_AC_h(X,param,XtX,Xty,Winit);
    [Yl,mu] = kmeans(X*w,k,'Replicates',param.nb_kms);
    Y = class2indicator(Yl,k);
    y_star = Y*mu;
    NormFrob(i+1) = trace((y_star-X*w)*(y_star-X*w)');
    Winit = w;
end
fprintf('\n');
Y = Yl;


%------- Display results :
nbGenes_fin = nb_Genes(w);
fprintf('number of genes before sparsity : %d \n',initSize);
fprintf('number of selected genes        : %d \n',nbGenes_fin);

XW = X*w;

% Silhouette values
[Si_sd,~] = silhouette(XW,Y);
fprintf('Silhouette average : %0.4f \n',mean(Si_sd));

% TSNE
if isTsne
    T=tsne(XW);
    figure('name','tsne : KSparse')
    gscatter(T(:,1),T(:,2),Y)
    title('tsne : KSparse');
end


end

