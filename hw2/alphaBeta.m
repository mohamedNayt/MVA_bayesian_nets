function [alpha, beta] = alphaBeta(init,trans,obs,means,covs)
    nx = size(trans,1);
    m = size(obs,1);
    %initialization
    alpha = zeros(nx,m);
    beta = zeros(nx,m);
    alpha(:,1) = log(obslik(obs(1,:),means,covs).*init);
    beta(:,m) = 0;
    %Recursion
    for t=2:m
        p=log(obslik(obs(t,:),means,covs));
        alpha(1,t) = logsum(p(1)+trans(:,1)+alpha(:,t-1));
        alpha(2,t) = logsum(p(2)+trans(:,2)+alpha(:,t-1));
        alpha(3,t) = logsum(p(3)+trans(:,3)+alpha(:,t-1));
        alpha(4,t) = logsum(p(4)+trans(:,4)+alpha(:,t-1));
        t2 = m-t+1;
        p=log(obslik(obs(t2+1,:),means,covs));
        beta(1,t2) = logsum(transpose(trans(1,:))+p+beta(:,t2+1));
        beta(2,t2) = logsum(transpose(trans(2,:))+p+beta(:,t2+1));
        beta(3,t2) = logsum(transpose(trans(3,:))+p+beta(:,t2+1));
        beta(4,t2) = logsum(transpose(trans(4,:))+p+beta(:,t2+1));
    end
end