function res=logp2(logalpha,logbeta,train,trans,means,covs)    
    %Given alpha beta compute (P(q(t)=j,q(t+1)=k|u(1)...u(T)))_{j,k,t}
    m=size(train,1);
    nx=size(trans,1);
    res = zeros(nx,nx,m-1);
    for t=1:m-1
        emis=repmat(log(transpose(obslik(train(t+1,:),means,covs))),[nx,1]);
        res(:,:,t) = log(trans)+emis+repmat(logalpha(:,t),[1,nx])+repmat(transpose(logbeta(:,t+1)),[nx,1]);
        temp = reshape(res(:,:,t),[nx*nx,1]);
        res(:,:,t) = res(:,:,t) - logsum(temp);
    end
end