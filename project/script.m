trans = [0.6,0.4;0.2,0.8];
obslik = [0.9,0.1;0.3,0.7];
init = [0.7;0.3];
m=10;
%%
[obs,states] = hmmgenerate(m,trans,obslik);
%%
[alpha,beta] = maxProduct(init,trans,obslik,obs);
%%
f = zeros(2,2,m-1);
for i=1:2
    for j=1:2
        for k=1:m-1
            f(i,j,k) = alpha(i,k)*beta(j,k+1)*obslik(i,j);
        end
    end
end

%% Finding the most likely sequence
x = zeros(3,m);
[mx,i] = max(f(:,:,1));
[~,j] = max(mx);
%%
x(1,1) = i(j);
x(1,2) = j;
for t=3:m
    [~,i] = max(f(x(1,t-1),:,t-1));
    x(1,t) = i;
end

%% Finding the second most likely sequence
n=2;
P = zeros(m,1);
out = f(:,:,1);
ind = 1:n;
out=out(ind~=x(1,1),:);
P(1) = max(out);

for i=2:m
    P(i) = max(f(x(1,i-1),ind~=x(1,i),i-1));
end
%%
[~,i0] = max(P);
if i0>1
    x(2,1:(i0-1)) = x(1,1:i0-1);
    out = f(x(2,i0-1),:,i0-1);
    out(x(1,i0)) = 0;
    [~,x(2,i0)] = max(out);
    for k=(i0+1):m
        [~,x(2,k)] = max(f(x(2,k-1),:,k-1));
    end
end
%% Finding the 3rd most likely sequence
P2 = zeros(m-i0+1,1);
P2(1)=0;
for i=2:(m-i0+1)
    out = f(x(2,i-1),:,i-1);
    out(x(2,i))=0;
    P2(i) = P(i0)*max(out)/f(x(2,i-1),x(2,i),i-1);
end

PP = [P,P2];
PP(i0)=0;
[~,i] = max(PP);
if i>m
    
else
    
end

%%
%% Test real data
obslik = xlsread('obslik.xls');
transit = xlsread('transit.xls');
transit(isnan(transit))=0;
obslik(isnan(obslik)) = 0;
%%
chars = ' aàâbcdeéèêëfghiîïjklmnoôpqrstuûvwxyz';
intchars = int16(chars);
char_2_int = zeros(2,37);
char_2_int(1,:) = intchars;
char_2_int(2,:) = 1:37;
%% normalize transitions matrix
transit = transit./repmat(sum(transit,2),[1 size(transit,2)]);
init = transpose(transit(1,:));
%% compute alpha and beta
[~,obs] = ismember(uint16('antecedent '),intchars);
[alpha,beta] = maxProduct(init,transit,obslik,obs);
% compute the f's
m=length(obs);
f = zeros(37,37,m-1);
for i=1:37
    for j=1:37
        for k=1:m-1
            f(i,j,k) = alpha(i,k)*beta(j,k+1)*transit(i,j)*obslik(j,obs(k+1));
        end
    end
end
% finding the most likely sequence
x = zeros(3,m);
[mx,i] = max(f(:,:,1));
[~,j] = max(mx);
x(1,1) = i(j);
x(1,2) = j;
for t=3:m
    [~,i] = max(f(x(1,t-1),:,t-1));
    x(1,t) = i;
end
%
decode = char(intchars(x(1,:)))
