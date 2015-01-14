function res=logsum(p)
    res=max(p)+log(sum(exp(p-repmat(max(p),[size(p,1),1]))));
end