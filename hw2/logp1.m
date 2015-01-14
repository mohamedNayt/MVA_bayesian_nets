function res=logp1(logalpha,logbeta)
%Given alpha beta compute (P(q(t)=k|u(1)...u(T)))_{t,k}
    res=logalpha+logbeta-repmat(logsum(logalpha+logbeta),[4,1]);
end