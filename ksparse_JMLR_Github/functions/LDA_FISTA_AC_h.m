function [W,NormL0,NormL1] = LDA_FISTA_AC_h(X,param,XtX,Xty,winit)
% Solves min norm(X*W - Y)^2 (frobenius norm) st norm(w,1) < param.LDA_ETA
% with Chambolle method

n = size(X,1);
d = size(X,2);

w_old = ones(d,param.h);
w_loc = winit;
t_old = 1;

NormL0 = zeros(1,param.LDA_MAXITER);
NormL1 = zeros(1,param.LDA_MAXITER);

for i = 1:param.LDA_MAXITER

    grad_w = XtX*w_loc - Xty;

    %gradient step
    w_new = w_loc - param.LDA_STEPSIZE*grad_w;

    w_new = reshape(w_new,d*param.h,1);

    %Projection on the l1 ball
    w_new_abs = abs(w_new)';
    x = zeros(size(w_new))';

    %mex fourni ball.zip
    proj_condat(w_new_abs,x,param.LDA_ETA);

    w_new = (sign(w_new)'.*x)';

    w_new((abs(w_new)<0.01),1)=0;

    %normL0 = sum(abs(w_new)>0.);

    NormL0(i) = sum(abs(w_new)>0);
    NormL1(i) = sum(abs(w_new));
    
    %Reshape back
    w_new = reshape(w_new,size(w_loc,1),size(w_loc,2));
    
    %Chambolle method
    t_new  = (i+5)/4;

    w_loc_new = w_new + ((t_old - 1)/t_new) * (w_new - w_old);

    w_old = w_new;
    w_loc = w_loc_new;
    t_old = t_new;

end

W = w_loc;

end
