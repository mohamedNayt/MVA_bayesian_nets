function [alpha, beta] = maxProduct(init,trans,obslik,obs)
    nx = size(trans,1);
    m = size(obs,2);

    alpha = zeros(nx,m);
    beta = zeros(nx,m);

    alpha(:,1) = obslik(:,obs(1)).*init;
    beta(:,m) = 1;

    for t=2:m
        alpha(:,t) = obslik(:,obs(t)).*transpose(max(trans.*repmat(alpha(:,t-1),1,nx)));
        t2 = m-t+1;
        beta(:,t2) = transpose(max(transpose(trans).*repmat(obslik(:,obs(t2+1)).*beta(:,t2+1),1,nx)));
    end
end