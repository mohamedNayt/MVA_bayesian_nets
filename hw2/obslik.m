function res=obslik(obs,means,covs)
    % Get P(q=k|obs)_{k}
    k=size(means,1);
    res=zeros(k,1);
    
    for j=1:k
        res(j)=mvnpdf(obs,means(j,:),covs(:,:,j));
    end
    
end