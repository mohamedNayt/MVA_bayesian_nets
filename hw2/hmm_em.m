function [init,trans,means,covs,logltr,loglts]=hmm_em(train,test,init,trans,means,covs,niter)
    
    eps=10^(-10);
    nx = size(trans,1);
    T = size(train,1);
    logltr=zeros(1,niter);
    loglts=zeros(1,niter);
    
    for j=1:niter
        %compute alpha beta and marginal probabilities on train and test
        %set: E step
        [logalpha,logbeta] = alphaBeta(init,trans,train,means,covs);
        [logalpha_ts,logbeta_ts] = alphaBeta(init,trans,test,means,covs);
        p1 = exp(logp1(logalpha,logbeta));
        p2 = exp(logp2(logalpha,logbeta,train,trans,means,covs));
        p1ts = exp(logp1(logalpha_ts,logbeta_ts));
        p2ts = exp(logp2(logalpha_ts,logbeta_ts,test,trans,means,covs));
        % M-Step
        init = p1(:,1);
        trans = sum(p2,3);
        trans = trans./repmat(sum(trans,2),[1,size(trans,2)]);
        means = p1*train;
        means = means./repmat(sum(p1,2),[1,2]);
        for k=1:nx
            temp = train-repmat(means(k,:),[T,1]);
            covs(:,:,k) = ((repmat(p1(k,:),[2,1]).*temp')*temp)./sum(p1(k,:));
        end
        % Compute likelihoods
        logltr(j)=init'*log(init+eps)+trace(sum(p2,3)*log(trans)');
        loglts(j)=init'*log(init+eps)+trace(sum(p2ts,3)*log(trans)');
        for k=1:nx
            temp = train-repmat(means(k,:),[T,1]);
            temp_ts = test-repmat(means(k,:),[T,1]);
            inv_cov = inv(covs(:,:,k));
            logltr(j)=logltr(j)-0.5*sum(sum((temp*inv_cov).*(temp.*repmat(p1(k,:)',1,2)))) - 0.5*sum(p2(k,:))*sum( log( eig( covs(:,:,k)) ) ) - log(2*pi);
            loglts(j)=loglts(j)-0.5*sum(sum((temp_ts*inv_cov).*(temp_ts.*repmat(p1ts(k,:)',1,2)))) - 0.5*sum(p2ts(k,:))*sum( log( eig( covs(:,:,k)) ) ) - log(2*pi);
        end
    end
end