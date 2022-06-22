%% ----------------------   SCRIPT CLUSTERING    --------------------------
fprintf('-----------------  SCRIPT CLUSTERING   ----------------------\n');
clear;
close all;
rng(10);

addpath(genpath('functions'));
addpath(genpath('Data_prepared'));

%-------------Import data :

%load('DataPatel.mat');
load('DataUsoskin.mat');
%load('DataKlein.mat');
%load('DataZeisel9.mat');


k = nb_clusters;

m = size(X,1);
d = size(X,2);


%-------------Ksparse :


%------KLEIN
%ETA_test = [150 5000 10000 15000 20000 30000 35000 46000];


%------ZEISEL
%ETA_test = [150 6000 12000 17000 20000 25000 30000 42000];
%ETA_test = [1000 6000 10000 12500 15000 17500 20000 22500 25000];
%ETA_test = [5000 9000 12000];


%------USOSKIN
ETA_test = [50 700 2000 4000 6000 8000 12000 16000 22000 28000];



% Advanced default parameters
param.LOOP = 10;               % number of loops
param.h = k+4;                 % Size of the centroids
param.nb_kms = 40;             % number of kmeans-replicates
param.sigma = 150;             % used for initialisation with spectral
initialization = 'spectral';   % 'spectral' or 'PCA'
param.LDA_MAXITER = 50;       % be careful 10 with FISTA 
param.LDA_STEPSIZE = 1/1.001;  % be careful 1/1.01 with Fista_AC 
isTsne = 0;                    % 1 to display tsne, 0 otherwise



SIL_ETA = zeros(1,size(ETA_test,2));
NBGENES_ETA = zeros(1,size(ETA_test,2));
ACC_ETA = zeros(1,size(ETA_test,2));
ARI_ETA = zeros(1,size(ETA_test,2));
NMI_ETA = zeros(1,size(ETA_test,2));

for f = 1:size(ETA_test,2)
    fprintf(['eta = ' num2str(ETA_test(f)) '\n']);
    param.LDA_ETA = ETA_test(f);         
    [Ysd,w,NormFrob,si_sd] = ksparse(X,k,param,initialization,isTsne);
    Ysd = match_names(YR,Ysd,k);
    SIL_ETA(f) = mean(si_sd);
    NBGENES_ETA(f) = nb_Genes(w);
    ACC_ETA(f) = compute_accuracy(YR,Ysd,k);
    ARI_ETA(f) = RandIndex(YR,Ysd);
    NMI_ETA(f) = nmi(YR,Ysd);
end




figure
subplot(2,1,1)
plot(ETA_test,SIL_ETA,'*-','LineWidth',2);
xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
xlim([ETA_test(1),ETA_test(end)])
ylim([min(SIL_ETA),1])
legend('Silhouette coefficient')
grid on
subplot(2,1,2)      
plot(NBGENES_ETA,SIL_ETA,'*-','LineWidth',2);
xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
xlim([NBGENES_ETA(1),NBGENES_ETA(end)])
ylim([min(SIL_ETA),1])
legend('Silhouette coefficient')
grid on






% figure
% subplot(2,1,1)
% plot(ETA_test,SIL_ETA,'*-','LineWidth',2);
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Silhouette coefficient')
% grid on
% subplot(2,1,2)      
% plot(ETA_test,ACC_ETA,'*-','LineWidth',2);
% hold on
% plot(ETA_test,ARI_ETA,'*-','LineWidth',2);
% hold on
% plot(ETA_test,NMI_ETA,'*-','LineWidth',2);
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on
% 
% 
% figure
% subplot(2,1,1)
% plot(NBGENES_ETA,SIL_ETA,'*-','LineWidth',2);
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
% legend('Silhouette coefficient')
% grid on
% subplot(2,1,2)      
% plot(NBGENES_ETA,ACC_ETA,'*-','LineWidth',2);
% hold on
% plot(NBGENES_ETA,ARI_ETA,'*-','LineWidth',2);
% hold on
% plot(NBGENES_ETA,NMI_ETA,'*-','LineWidth',2);
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on




% figure
% subplot(2,1,1)
% xx = ETA_test(1):20:ETA_test(end);
% yySil = interp1(ETA_test,SIL_ETA,xx,'pchip');
% plot(xx,yySil,'LineWidth',2)
% xlabel('constraint \eta')
% legend('Silhouette coefficient')
% grid on
% subplot(2,1,2)
% yyAcc = interp1(ETA_test,ACC_ETA,xx,'pchip');
% plot(xx,yyAcc,'LineWidth',2)
% hold on
% yyAri = interp1(ETA_test,ARI_ETA,xx,'pchip');
% plot(xx,yyAri,'LineWidth',2)
% hold on
% yyNmi = interp1(ETA_test,NMI_ETA,xx,'pchip');
% plot(xx,yyNmi,'LineWidth',2)
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on
% 
% 
% 
% figure
% subplot(2,1,1)
% xx = NBGENES_ETA(1):20:NBGENES_ETA(end);
% yySil = interp1(NBGENES_ETA,SIL_ETA,xx,'pchip');
% plot(xx,yySil,'LineWidth',2)
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold')
% legend('Silhouette coefficient')
% grid on
% subplot(2,1,2)
% yyAcc = interp1(NBGENES_ETA,ACC_ETA,xx,'pchip');
% plot(xx,yyAcc,'LineWidth',2)
% hold on
% yyAri = interp1(NBGENES_ETA,ARI_ETA,xx,'pchip');
% plot(xx,yyAri,'LineWidth',2)
% hold on
% yyNmi = interp1(NBGENES_ETA,NMI_ETA,xx,'pchip');
% plot(xx,yyNmi,'LineWidth',2)
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on
% 



% figure
% [NBGENES_ETA_SORT, indSort] = sort(NBGENES_ETA);
% subplot(2,1,1)
% plot(NBGENES_ETA_SORT,SIL_ETA(indSort),'*-','LineWidth',2);
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
% legend('Silhouette coefficient')
% grid on
% subplot(2,1,2)      
% plot(NBGENES_ETA_SORT,ACC_ETA(indSort),'*-','LineWidth',2);
% hold on
% plot(NBGENES_ETA_SORT,ARI_ETA(indSort),'*-','LineWidth',2);
% hold on
% plot(NBGENES_ETA_SORT,NMI_ETA(indSort),'*-','LineWidth',2);
% xlabel('Number of selected genes','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on



% figure
% subplot(3,1,1)
% plot(ETA_test,SIL_ETA,'*-','LineWidth',2);
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Silhouette coefficient')
% grid on
% subplot(3,1,2)      
% plot(ETA_test,ACC_ETA,'*-','LineWidth',2);
% hold on
% plot(ETA_test,ARI_ETA,'*-','LineWidth',2);
% hold on
% plot(ETA_test,NMI_ETA,'*-','LineWidth',2);
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Accuracy','ARI','NMI')
% grid on
% subplot(3,1,3)
% plot(ETA_test,NBGENES_ETA,'*-','LineWidth',2);
% xlabel('constraint \eta','FontSize',12,'FontWeight','bold');
% legend('Number of selected genes')
% grid on






