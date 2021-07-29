function [w, d, S, Q] = scTRCA(X, Xb)
% The computation process refers to literature [1]
% [1]	H. Tanaka, ¡°Group task-related component analysis (gTRCA): 
%       a multivariate method for inter-trial reproducibility and inter-subject 
%       similarity maximization for EEG data analysis,¡± Sci Rep, vol. 10, no. 1,
%       pp. 84, Jan 9, 2020.
% Inputs:
%       X: data in cell format, X = cell(1, num)
%       X{num} : Nchannels x Nsamples
%       Xb{num}: Nchannels x tau x Ntrials

[~, Num] = size(X);
[~, tau, ~] = size(Xb{1});

%% computation of U, V, Q0 matrices:
U = cell(Num,1);
V = cell(Num,1);
Q0 = cell(Num,1);
for n=1:Num
    [Nc, Np, ~] = size(Xb{n});
    U{n}=zeros(Nc, Np);
    V{n}=zeros(Nc, Nc);
    Q0{n}=zeros(Nc, Nc);
end

for n = 1:Num
    % U
    U{n} = mean(Xb{n},3);
    % V
    K = size(Xb{n},3);
    for k=1:K
        V{n} = V{n} + Xb{n}(:,:,k)*Xb{n}(:,:,k)';
    end
    V{n} = V{n}/K;
    % Q
    [~, T] = size(X{n});
    Q0{n} = X{n}*X{n}'/T;
end

%% computation of S and Q matrices:
% S matrix
S = cell(2);
K1 = size(Xb{1},3);
S{1,1} = K1/((K1-1)*tau)*(U{1}*U{1}' - V{1}/K1);
S{1,2} = 1/tau*(U{1}*U{2}');
S{2,1} = 1/tau*(U{2}*U{1}');
S{2,2} = 1/tau*(U{2}*U{2}');
% S{2,2} = zeros(size(S{2,1},1),size(S{1,2},2));
S = cell2mat(S);

% Q matrix
Q = blkdiag(Q0{1}, Q0{2});

%% generalized eigendecomposition:
[W, D] = eig(Q\S); 
D = diag(D);
[~, index] = sort(D, 'descend');
d = D(index(1));
w = W(:,index(1));


