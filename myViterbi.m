function path = myViterbi(transMat, loglikeMat, initProb)
% Implementation of the Viterbi algorithm to find the path of states that has the
% highest posterior probability in a finite-state hidden Markov model
% (HMM).
%
% Input
%   - transMat      : the transition probability matrix, P(S_n | S_n-1).
%                       Rows correspond to the starting states; columns
%                       corresond to the ending states. Each row sums to 1.
%                       (size: nState * nState)
%   - loglikeMat    : the log-likelihood matrix of each state to explain
%                       each observation, log( P(O_n | S_n) ) Rows
%                       correspond to states; columns correspond to
%                       observations. (size: nState * nObserv)
%   - initProb      : the initial probability of all states, P(S_0). A
%                       column vector with nState elements.
%
% Output
%   - path          : the best path of states (S_1, ..., S_N) that gives
%                       the maximal posterior probability, P(S_1, ..., S_N
%                       | O_1, ... O_N). A column vector with nObserv elements.
%


if nargin<3     % use a uniform distribution if the initial probability is not given
    initProb = ones(size(transMat,1),1)/size(transMat,1); 
end
if size(transMat,1) ~= size(loglikeMat, 1) || size(transMat,1) ~= length(initProb)
    error('The number of states is not consistent in the transition matrix, the likelihood matrix, and the initial probability!');
end

% Your implementation starts here...
N = size(transMat,1); % state space size
T = size(loglikeMat,2); % sequence length 
S = 1:1:size(loglikeMat,1);

v(1,:) = log(initProb)+loglikeMat(:,1); % initialize v_1(j)
prev(1,:) = 0;
for t = 2:T
    for j = 1:N
        v(t,j) = max(v(t-1,:)+log(transMat(:,j)'))+loglikeMat(j,t);
        [~,prev(t,j)] = max(v(t-1,:)+log(transMat(:,j)'));
    end
end
[~,idx(T)] = max(v(T,:)); 
path(T) = S(idx(T));
for t = T:-1:2
    idx(t-1) = prev(t,idx(t));
    path(t-1) = S(idx(t-1));
end
end

