
 function [w,mu]=FISTA_Primal(X,YR,k,param,isEpsilon)
% Solves min norm(X*W - Y)^2 (frobenius norm) st norm(w,1) < param.eta
% with Chambolle method

n = size(X,1);
d = size(X,2);
Y = class2indicator(YR,k);
XtX=X'*X;
Xty=X'*Y;

w_old = ones(d,k);
w_loc = w_old;
t_old = 1;


for i = 1:param.niter

    grad_w = XtX*w_loc - Xty;

    %gradient step
    V = w_loc - param.gamma*grad_w;

    V= reshape(V,d*k,1);

    %Projection on the l1 ball
    proj_l1ball = @(y) max(abs(y)-max(max((cumsum(sort(abs(y),1,'descend'),1)-param.eta)./(1:size(y,1))'),0),0).*sign(y);
    
       V=proj_l1ball(V);
    
%        V((abs(V)<0.01),1)=0;
    
    %Reshape back
   
        w_new = reshape(V,d,k);
    
   
    %Chambolle method
    t_new  = (i+5)/4;

    w_loc_new = w_new + ((t_old - 1)/t_new) * (w_new - w_old);

    w_old = w_new;
    w_loc = w_loc_new;
    t_old = t_new;

end

w = w_loc;
mu = centroids(X*w,YR,k);

   nbGenes_fin= nb_Genes(w);

end
