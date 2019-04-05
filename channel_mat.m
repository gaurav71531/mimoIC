function channel_mat(n_t, n_r, users)
% The code is written to reproduce the results of paper:

% Gaurav Gupta and A.K. Chaturvedi, "User Selection in MIMO Interfering
% Broadcast Channels", IEEE TCOM 2014

% The code is open-access, and for using the code the user is required to 
% cite the above paper.


% n_t = 4;
% n_r = 6;
% users = 2;
N_sim = 1000;
directive = {'MU_MIMO', 'MIMO_IC'};
command = 2;
BS = 2;

chDir = 'Channel_Mat/';

if ~exist(chDir, 'dir')
    mkdir(chDir);
end
    

for i = users
    switch command
        case 1
            H_all = (randn(n_r, n_t,i,N_sim) + sqrt(-1) * randn(n_r, n_t,i,N_sim)) / sqrt(2);
            str = sprintf('H_all_%d_%d_%d',n_t,n_r,i);
        case 2
            H_all = (randn(n_r, n_t,i,BS,BS,N_sim) + sqrt(-1) * randn(n_r, n_t,i,BS,BS,N_sim)) / sqrt(2);
            str = sprintf('H_all_%d_%d_%d_%d',n_t,n_r,i,BS);
    end
    save([chDir str], 'H_all');
end