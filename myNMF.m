function [W, H, KL] = myNMF(V, r, nIter, bUpdateW, bUpdateH,initW, initH)
% Implementation of the multiplicative update rule for NMF with K-L
% divergence. W should always be normalized so that each column of W sums
% to 1. Scale H accordingly. This normalization step should be performed in
% each iteration, and in initialization as well. This is to deal with the
% non-uniqueness of NMF solutions.
%
% Input
%   - V         : a non-negative matrix (m*n) to be decomposed
%   - r         : #columns of W (i.e., #rows of H)
%   - nIter     : #iterations
%   - initW     : initial value of W (m*r) (default: a random non-negative matrix)
%   - initH     : initial value of H (r*n) (default: a random non-negative matrix)
%   - bUpdateW  : update W (bUpdateW==1) or not (bUpdateW==0) (default: 1)
%   - bUpdateH  : update H (bUpdateH==1) or not (bUpdateH==0) (default: 1)
%
% Output
%   - W         : learned W
%   - H         : learned H
%   - KL        : KL divergence after each iteration (a row vector with nIter+1
%               elements). KL(1) is the KL divergence using the initial W
%               and H.
%
% Author: Zhiyao Duan
% Created: 10/9/2013
% Last modified: 10/2/2014

[m,n] = size(V);
if nargin<7 initH = rand(r, n); end     % randomly initialize H
if nargin<6 initW = rand(m, r); end     % randomly initialize W
if nargin<5 bUpdateH=1; end
if nargin<4 bUpdateW=1; end
if r~=size(initW,2) || r~=size(initH,1)
    error('Parameter r and the size of W or H do not match!');
end

% Your implementation starts here...
% step 1: initialization of W and H according to initW & initH
% First, check if there is any zero element in initW and initH
if nargin<6 % No inital W input -> when randomly intialize W
    while (initW == 0) % outputs 1 if there is zero in initW 
        initW = rand(m,r); % Makes sure there is no zero element
    end
end
if nargin<7 % No inital H input -> when randomly intialize H
    while (initH == 0) % outputs 1 if there is zero in initH
        initH = rand(r,n); % Makes sure there is no zero element
    end
end
% solve non-uniquness issue by normalization.
W = initW./sum(initW,1); 
H = initH./sum(initW,1)';
I = ones(m,n);
% step 2: start a loop 
Pd = W*H; % Reconstruct V for the first element
KL(1) = sum(V(:).*log(V(:)./Pd(:)) - V(:) + Pd(:)); % First KL Div.
for n = 1:nIter
    if bUpdateW == 1 % Update W based on multiplicative rule (KL Div.)
        W = W.*((V./Pd)*H')./(I*H');
        W = W./sum(W,1); % Normalization
    end
    Pd = W*H; % Reconstruct V after updating W
    if bUpdateH == 1 % Update H based on multiplicative rule (KL Div.)
        H = H.*((W'*(V./Pd))./(W'*I));
        if bUpdateW == 1 % When W gets updated we need to normalize H as well
            H = H./sum(W,1)'; % Normalization
        end
    end
   Pd = W*H; % Reconstruct V after updating H
   KL(n+1) = sum(V(:).*log(V(:)./Pd(:)) - V(:) + Pd(:)); % Update KL Div.
end

end