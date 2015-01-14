function [seq,p_max] = viterbi(init,trans,obs,means,covs)
    nx = size(trans,1);
    m = size(obs,1);
    seq=zeros(1,m);
    arg_alpha = zeros(nx,m);
    alpha = zeros(nx,m);
    alpha(:,1) = log(obslik(obs(1,:),means,covs).*init);
    %Compute max-product alpha and the transition led to the best sequence
    %at t
    for t=2:m
        p=log(obslik(obs(t,:),means,covs));
        alpha(1,t) = max(p(1)+trans(:,1)+alpha(:,t-1));
        alpha(2,t) = max(p(2)+trans(:,2)+alpha(:,t-1));
        alpha(3,t) = max(p(3)+trans(:,3)+alpha(:,t-1));
        alpha(4,t) = max(p(4)+trans(:,4)+alpha(:,t-1));
        [~,arg_alpha(1,t)] = max(trans(:,1)+alpha(:,t-1));
        [~,arg_alpha(2,t)] = max(trans(:,2)+alpha(:,t-1));
        [~,arg_alpha(3,t)] = max(trans(:,3)+alpha(:,t-1));
        [~,arg_alpha(4,t)] = max(trans(:,4)+alpha(:,t-1));
    end
    
    [p_max,seq(m)]=max(alpha(:,m));
    %Decode
    for t=(m-1):-1:1
        seq(t)=arg_alpha(seq(t+1),t+1);
    end
    
end