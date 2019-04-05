clear all

% The code is written to reproduce the results of paper:

% Gaurav Gupta and A.K. Chaturvedi, "User Selection in MIMO Interfering
% Broadcast Channels", IEEE TCOM 2014

% The code is open-access, and for using the code the user is required to 
% cite the above paper.
% 
%% generate channel matrices

n_t = 3; % number of transmit antennas
n_r = 2; % number of receiver antennas
users = 5:5:40;  % x-axis of plot (total number of users in the system)
BS = 2;

channel_mat(n_t, n_r, users);
SNR = 20; % in dB

% execute algorithms

%% 1. c-algorithm

MIMO_IC('group_cap', n_t, n_r, SNR, BS, users);
% use dof = 1
%% 2. o-algorithm

MIMO_IC('group_cho', n_t, n_r, SNR, BS, users);
% use dof = 1
%% 2. single-user

MIMO_IC('ic_cap', n_t, n_r, SNR, BS, users);

%% plot results

figure;

plot(users, R_F_grp_cap, users, R_F_grp_cho1, users, R_F_IC_cap);
grid;
