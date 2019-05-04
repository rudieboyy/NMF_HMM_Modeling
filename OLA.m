function y = OLA(X,win,M,R)
% Input
% X - spectrogram (time*freq) matrix
% win - window function
% M - length of a window
% R - hop size

% Output 
% y - reconstructed output 

bins = size(X,1);
nframes = size(X,2); % number of frames
X = X(:);
time = nframes*R+M; %time = nframes (loops) * hop size + window length
y = zeros(time,1);
n = 1;
for m = 1:nframes
    X_temp = X(m*bins-(bins-1):m*bins); % take one frame
    X_temp = [X_temp;conj(X_temp(end-1:-1:2))]; %concatenate with conjugate values (negative freq comp) 
    y_temp = real(ifft(X_temp)); % take ifft of one frame
    y_temp = y_temp(1:M); % take window length y[n']
    y(n:(n+M-1)) = y(n:(n+M-1))+(y_temp.*win); %sum y and y[n']*w[n'-mR]
    n = n+R; % move the index by hop size
end

end