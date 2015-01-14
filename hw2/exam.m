train=importdata('EMGaussian.data');
test=importdata('EMGaussian.test');
%% Values from prrevious homework
means=[-3.06196072, -3.53454048;-2.03436695,  4.17258596; 3.97793025,  3.77333059;3.80070949, -3.79729742];
a= [6.24140909,  6.05017464; 6.05017464,  6.18245528];
b= [2.90442381, 0.20655763; 0.20655763,  2.75617077];
c= [0.21035667, 0.29045085; 0.29045085,  12.23996609];
d= [0.92127927, 0.0573808; 0.0573808, 1.86586017];
covs=cat(3,a,b,c,d);
init= ones(4,1)/4;
A=ones(4,4)/6+diag(ones(1,4)/3);
%% Compute the alpha beta function by recursion on the training set
[logalpha,logbeta]=alphaBeta(init,log(A),train,means,covs);
%% We compute and plot the marginal probabilities using alpha and beta
lp1 = logp1(logalpha,logbeta);
figure
subplot(2,2,1)
plot(exp(lp1(1,1:100)),'o')
title('P(q(t)=1|u(1)...u(T))');
subplot(2,2,2)
plot(exp(lp1(2,1:100)),'o');
title('P(q(t)=2|u(1)...u(T))');
subplot(2,2,3)
plot(exp(lp1(3,1:100)),'o');
title('P(q(t)=3|u(1)...u(T))');
subplot(2,2,4)
plot(exp(lp1(4,1:100)),'o');
title('P(q(t)=4|u(1)...u(T))');
%% We compute the marginal probability P(q(t),q(t+1)|u(1)...u(T)) represented as matrix
lp2 = logp2(logalpha,logbeta,train,A,means,covs);
%% HMM EM performs the EM algorithm
[init2,trans2,means2,covs2,logltr,loglts]=hmm_em(train,test,init,A,means,covs,10);
%% Plot the log likelihoods
figure
subplot(1,2,1)
plot(logltr);
title('loglikelihood vs time for training data');
subplot(1,2,2)
plot(loglts);
title('loglikelihood vs time for test data');
%% Plot training data and clusters centers
scatter(train(:,1),train(:,2));hold on,
scatter(means2(:,1),means2(:,2),'rx','LineWidth',2);
%% The function viterbi implements the viterbi algorithm
[seq,p_max]=viterbi(init2,trans2,train,means2,covs2);
%% Training data + clusters centers + Viterbi classification
scatter(train(seq==1,1),train(seq==1,2),'ro');hold on,
scatter(train(seq==2,1),train(seq==2,2),'bo');hold on,
scatter(train(seq==3,1),train(seq==3,2),'go');hold on,
scatter(train(seq==4,1),train(seq==4,2),'yo');hold on,
scatter(means2(:,1),means2(:,2),'kx','LineWidth',2);
%%
scatter(test(:,1),test(:,2));hold on,
scatter(means2(:,1),means2(:,2),'rx','LineWidth',2);
%% Compute alpha beta and marginal probabilities for the test set
[logalpha,logbeta]=alphaBeta(init2,log(trans2),test,means2,covs2);
lp1ts = logp1(logalpha,logbeta);
figure
subplot(2,2,1)
plot(exp(lp1ts(1,1:100)),'o')
title('P(q(t)=1|u(1)...u(T))');
subplot(2,2,2)
plot(exp(lp1ts(2,1:100)),'o')
title('P(q(t)=2|u(1)...u(T))');
subplot(2,2,3)
plot(exp(lp1ts(3,1:100)),'o')
title('P(q(t)=3|u(1)...u(T))');
subplot(2,2,4)
plot(exp(lp1ts(4,1:100)),'o')
title('P(q(t)=4|u(1)...u(T))');
%% Compute the most likely state for each marginal probability+plot
[~,likely_states] = max(lp1ts);
plot(likely_states(1:100),'ro');
title('Likely state vs time');
%%

%% Implement the Viterbi algorithm on test data and compare it with the previous approach
[seq,p_max]=viterbi(init2,trans2,test,means2,covs2);
plot(seq(1:100),'ro');hold on,
plot(likely_states(1:100),'b+');
legend('Viterbi decoding','Most likely state');
title('States inference using Viterbi and maximum marginal probability');
%% Plot test set + clusters centers + most likely sequence
scatter(test(seq==1,1),test(seq==1,2),'ro');hold on,
scatter(test(seq==2,1),test(seq==2,2),'bo');hold on,
scatter(test(seq==3,1),test(seq==3,2),'go');hold on,
scatter(test(seq==4,1),test(seq==4,2),'yo');hold on,
scatter(means2(:,1),means2(:,2),'kx','LineWidth',2);
%% Compare most ikely sequence and most likely state for each marginal probability
scatter(test(and(likely_states==1,likely_states==seq),1),test(and(likely_states==1,likely_states==seq),2),'ro');hold on,
scatter(test(and(likely_states==2,likely_states==seq),1),test(and(likely_states==2,likely_states==seq),2),'bo');hold on,
scatter(test(and(likely_states==3,likely_states==seq),1),test(and(likely_states==3,likely_states==seq),2),'go');hold on,
scatter(test(and(likely_states==4,likely_states==seq),1),test(and(likely_states==4,likely_states==seq),2),'yo');hold on,

scatter(test(and(likely_states==1,likely_states~=seq),1),test(and(likely_states==1,likely_states~=seq),2),'ro','LineWidth',3);hold on,
scatter(test(and(likely_states==2,likely_states~=seq),1),test(and(likely_states==2,likely_states~=seq),2),'bo','LineWidth',3);hold on,
scatter(test(and(likely_states==3,likely_states~=seq),1),test(and(likely_states==3,likely_states~=seq),2),'go','LineWidth',3);hold on,
scatter(test(and(likely_states==4,likely_states~=seq),1),test(and(likely_states==4,likely_states~=seq),2),'yo','LineWidth',3);hold on,

scatter(means2(:,1),means2(:,2),'kx','LineWidth',2);