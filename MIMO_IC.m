function MIMO_IC(varargin)

% % % 
% The code is written to reproduce the results of paper:

% Gaurav Gupta and A.K. Chaturvedi, "User Selection in MIMO Interfering
% Broadcast Channels", IEEE TCOM 2014

% The code is open-access, and for using the code the user is required to 
% cite the above paper.
% % % 
% The following package is used for visual ease:
% progressbar, % Author:Steve Hoelzer
% % % 

% sample input : (directive, n_t, n_r, SNR, BS, users)
% directive : code of algorithm to implement (see details of code below)
% n_t, n_r, SNR as usual
% BS : no of base stations comprising MIMO Interference channel
% users : users number vector
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

global n_t
global n_r
ii = 2;
n_t = varargin{ii};
n_r = varargin{ii+1};
SNR = varargin{ii+2};
BS = varargin{ii+3};
users_num = varargin{ii+4};

switch lower(varargin{1})
    
    case 'ic'
        R_F = SumCapacity_MIMO_IC_single_user(SNR, BS, users_num);
        str = 'R_F_IC';
        
    case 'ic_cap'
        R_F = SumCapacity_MIMO_IC_capacity(SNR, BS, users_num);
        str = 'R_F_IC_cap';
    
    case 'optimal'
        R_F = SumCapacity_optimal(SNR, BS, users_num);
        str = 'R_F_opt_IC';
        
    case 'capacity'
        R_F = SumCapacity_capacity(SNR, BS, users_num);
        str = 'R_F_cap_IC';
        
    case 'conde'
        R_F = SumCapacity_conditional_Ent(SNR, BS, users_num);
        str = 'R_F_gg_IC';
        
    case 'convention'
        R_F = SumCapacity_conventional(SNR, BS, users_num);
        str = 'R_F_conv_IC';
        
    case 'tdma'
        R_F = SumCapacity_TDMA(SNR, BS, users_num);
        str = 'R_F_TDMA_IC';
        
    case 'dist_iter'
%         Here number of users are number of transmitter,Receiver pairs
        R_F = SumCapacity_distributed_iterative(SNR, BS);
        str = 'R_F_Diter';
        
    case 'max_sinr'
        R_F = SumCapacity_iterative_max_SINR(SNR, BS);
        str = 'R_F_mSINR';
        
    case 'group'
        R_F = SumCapacity_grouping(SNR, BS, users_num);
        str = 'R_F_grp';
        
    case 'group_cap'
        R_F = SumCapacity_grouping_Capacity(SNR, BS, users_num);
        str = 'R_F_grp_cap';
        
    case 'cap_test'
        R_F = SumCapacity_grouping_CapacityTest(SNR, BS, users_num);
        str = 'R_F_grp_capT';
        
    case 'group_cho'
        R_F = SumCapacity_grouping_Chordal(SNR, BS, users_num);
        str = 'R_F_grp_cho1';
        
    case 'group_opt'
        R_F = SumCapacity_grouping_Optimal(SNR, BS, users_num);
        str = 'R_F_grp_opt';
        
    case 'three_cell'
        R_F = SumCapacity_three_cell_alignment(SNR, BS, users_num);
        str = 'R_F_3_cell';
        
    case 'three_cell_cap'
        R_F = SumCapacity_three_cell_Capacity(SNR, BS, users_num);
        str = 'R_F_3_cell_cap';
        
    case 'three_cell_cho'
        R_F = SumCapacity_three_cell_Chordal(SNR, BS, users_num);
        str = 'R_F_3_cell_cho';
        
    case 'm2'
%         select different models in the function itself
        [R_F,appr] = SumCapacity_cyclic_model2(SNR, BS, users_num);
        str = ['R_F_M2_' appr];
        
    case 'm2_cap'
%         select different models in the function itself
        [R_F,appr] = SumCapacity_cyclic_model2_Capacity(SNR, BS, users_num);
        str = ['R_F_M2_' appr '_cap'];
        
    case 'hybrid'
        R_F = SumCapacity_hybrid_grouping(SNR, BS, users_num);
        str = 'R_F_hyb';
        
    case 'gg'
        R_F = SumCapacity_GG(SNR, BS, users_num);
        str = 'R_F_gg';
        
    case 'test'
        R_F = SumCapacity_Test(SNR, BS, users_num);
        str = 'R_F_test';
        
    otherwise
        disp('Wrong directive input, try again!');
        return;
end
assignin('base', str, R_F);


function R_F = SumCapacity_MIMO_IC_single_user(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
N_sim = 1000;
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        comb = ones(2,1);
        for avg_sim_ind = 1:N_sim

%             H = (randn(n_r,n_t,2,2,nusers) + 1i*randn(n_r,n_t,2,2,nusers))/sqrt(2);
            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>         
            R_final(avg_sim_ind) = Sum_capacity_MIMO_IC(H, SNR_vect(SNR_ind),comb);
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end 
        R_F(SNR_ind, mi) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_conditional_Ent(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            K_set = zeros(1,BS);
            for cell_ind = 1:BS
                Frob_temp = 0;
                for user_ind = 1:nusers
                    Frob_Norm = norm(H(:,:,user_ind,cell_ind,cell_ind), 'fro')^2;
                    if Frob_Norm > Frob_temp
                        Frob_temp = Frob_Norm;
                        maxValInd = user_ind;
                    end
                end
                K_set(cell_ind) = maxValInd;
            end
            SNR_tmp = 10^(SNR_vect(SNR_ind)/10) / n_t;
            for i = 1:BS
                Ent_temp = 0;
                for j = 1:nusers
                    H1 = H(:,:,j,i,i);
                    H2 = [];
                    for ui = 1:BS
                        if ui == i, continue;end
                        if ~isempty(H2)
                            H3 = H(:,:,j,i,ui);
                            break;
                        end
                        H2 = H(:,:,j,i,ui);
                    end
                    Ent123 = log2(real(det(eye(n_t) + SNR_tmp*((H1'*H1) + (H2'*H2) + (H3'*H3)))));
                    Ent23 = log2(real(det(eye(n_t) + SNR_tmp*((H2'*H2) + (H3'*H3)))));
                    Ent = Ent123 - Ent23;
                    T = eye(n_t);
%                     for k = 1:3
%                         if k == i, continue;end
%                         H1 = H(:,:,K_set(k), k, k);
%                         H2 = [];
%                         for ui = 1:BS
%                             if ui == k, continue;end
%                             if ~isempty(H2),
% %                                 T = H(:,:,j,i,k) \ H(:,:,j,i,ui);
% %                                 T = H(:,:,j,i,ui) \ H(:,:,j,i,k);
% %                                 T = inv(H(:,:,j,i,ui));
% %                                 T = inv(H(:,:,j,ui,k));
% %                                 ind_3 = ui;
%                                 H3 = T * H(:,:,K_set(k),k,ui);
% %                                 H3 = H(:,:,K_set(k),k,ui);
%                                 break;
%                             end
% %                             T = inv(H(:,:,j,i,ui));
% %                             T = inv(H(:,:,j,ui,k));
% %                             ind_2 = ui;
%                             H2 = T * H(:,:,K_set(k),k,ui);
%                         end 
% %                         T = H(:,:,j,i,ind_2) \ H(:,:,j,i,ind_3);
% %                         H2 = T * H2;
% %                         T = H(:,:,j,i,ind_3) \ H(:,:,j,i,ind_2);
% %                         H3 = T*H3;
% %                         H1 = T * H1;
% %                         H2 = T * H2;
%                         Ent123 = log2(real(det(eye(n_t) + SNR_tmp*((H1'*H1) + (H2'*H2) + (H3'*H3)))));
%                         Ent23 = log2(real(det(eye(n_t) + SNR_tmp*((H2'*H2) + (H3'*H3)))));
%                         Ent = Ent + Ent123 - Ent23;
% %                         Ent = Ent + Ent123;
%                     end
%                     sz = 1;
%                     mat = cell(1,2);
%                     for cell_ind = 1:BS
%                         if cell_ind == i, continue;end
%                         mat{sz} = H(:,:,j,i,cell_ind)'*H(:,:,j,i,cell_ind)*SNR_tmp;
%                         sz = sz + 1;
%                     end
%                     Ent12 = log2(real(det(eye(n_t) + mat{1} + mat{2})));
%                     Ent1 = log2(real(det(eye(n_t) + mat{1})));
%                     Ent2 = log2(real(det(eye(n_t) + mat{2})));
%                     Ent = 2*Ent12 - Ent1 - Ent2;
%                     Ent = log2(real(det(eye(n_t) + H(:,:,j,i,i)'*H(:,:,j,i,i)*SNR_tmp)));
%                     K_set_tmp = K_set;
%                     K_set_tmp(i) = j;
%                     [U, V] = compute_post_pre_processing_matrices_MMSE_opt(H, K_set_tmp, SNR_vect(SNR_ind));
%                     A = [];
%                     for ui = 1:BS
%                         if ui == i, continue;end
%                         if ~isempty(A),
%                             B = H(:,:,j,i,ui)*V{ui};
%                             break;
%                         end
%                         A = H(:,:,j,i,ui)*V{ui};
%                     end
%                     Ent = cho_dist(A, B);
                    if Ent > Ent_temp
                        Ent_temp = Ent;
                        maxValInd = j;
                    end
                end
                K_set(i) = maxValInd;
            end
            R_final(avg_sim_ind) = sum_capacity_R_sigma(H, SNR_vect(SNR_ind), K_set);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', users_num(mi), SNR_vect(SNR_ind), time(SNR_ind,mi));
    end 
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_capacity(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            K_set = zeros(1,BS);
            for cell_ind = 1:BS
                Frob_temp = 0;
                for user_ind = 1:nusers
                    Frob_Norm = norm(H(:,:,user_ind,cell_ind,cell_ind), 'fro')^2;
                    if Frob_Norm > Frob_temp
                        Frob_temp = Frob_Norm;
                        maxValInd = user_ind;
                    end
                end
                K_set(cell_ind) = maxValInd;
            end

            for cell_ind = 1:BS
                R_temp = 0;
                for user_ind = 1:nusers
                    K_set_temp = K_set;
                    K_set_temp(cell_ind) = user_ind;
                    R_S = sum_capacity_R_sigma(H, SNR_vect(SNR_ind), K_set_temp);
                    if R_S > R_temp
                        R_temp = R_S;
                        maxValInd = user_ind;
                    end
                end
                K_set(cell_ind) = maxValInd;
                R_final(avg_sim_ind) = max(R_final(avg_sim_ind),R_temp);
            end
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', users_num(mi), SNR_vect(SNR_ind), time(SNR_ind,mi));
    end 
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_MIMO_IC_capacity(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
N_sim = 1000;
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        for avg_sim_ind = 1:N_sim

%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            H = (randn(n_r,n_t,2,2,nusers) + 1i*randn(n_r,n_t,2,2,nusers))/sqrt(2);
            K_set = zeros(1,2);
            for cell_ind = 1:2
                Frob_temp = 0;
                for user_ind = 1:nusers
                    Frob_Norm = norm(H(:,:,cell_ind,cell_ind,user_ind), 'fro')^2;
                    if Frob_Norm > Frob_temp
                        Frob_temp = Frob_Norm;
                        maxValInd = user_ind;
                    end
                end
                K_set(cell_ind) = maxValInd;
            end
            
            R_temp = 0;
            maxValInd_exist = 0;
            for cell_ind = 1:2
                for user_ind = 1:nusers
                    K_set_temp = K_set;
                    K_set_temp(cell_ind) = user_ind;
                    R_S = Sum_capacity_MIMO_IC(H, SNR_vect(SNR_ind), K_set_temp);
                    if R_S > R_temp
                        R_temp = R_S;
                        maxValInd = user_ind;
                        maxValInd_exist = 1;
                    end
                end
                if maxValInd_exist
                    K_set(cell_ind) = maxValInd;
                end
                maxValInd_exist = 0;
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', users_num(mi), SNR_vect(SNR_ind), time(SNR_ind,mi));
    end 
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_optimal(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            R_temp = 0;
            for i = 1:nusers
                for j = 1:nusers
                    for k = 1:nusers
                        R_S = sum_capacity_R_sigma(H, SNR_vect(SNR_ind), [i,j,k]);
                        if R_S > R_temp
                            R_temp = R_S;
                        end
                    end
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_hybrid_grouping(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            R_final(avg_sim_ind) = Sum_capacity_hybrid_grouping(H, SNR_vect(SNR_ind), nusers, BS, DoF);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_S = Sum_capacity_hybrid_grouping(H, SNR, K, L, dof)

global n_t
global n_r

SNR = 10^(SNR/10);
sz = n_t - ((K(L-2)+1)*dof);
U = zeros(n_r,dof,K,L);
V = zeros(n_t,sz,K,L);
G = zeros(n_t,dof,L);

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if cell_ind + 1 > L, next_cell_ind = 1;end
    if K <= 3
        mat_req = zeros(K*n_t,n_t + K*n_r);
        sz_eye = 1;
        sz_H_colm = n_t;
        sz_H_row = 1;
        for i = 1:K
            mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
            mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,i,next_cell_ind,cell_ind)';
            sz_eye = sz_eye + n_t;
            sz_H_row = sz_H_row + n_t;
            sz_H_colm = sz_H_colm + n_r;
        end
        X = null(mat_req);
        G(:,:,cell_ind) = X(1:n_t,:);
        sz = n_t;
        for i = 1:K
            U(:,:,i,next_cell_ind) = X(sz+1:sz+n_r,:);
            sz = sz + n_r;
        end
    else
        disp('bye');
    end
end

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if next_cell_ind > L, next_cell_ind = 1;end
    for user_ind = 1:K
        mat_req = zeros((K*(L-2)+1)*dof,n_t);
        mat_req(1:dof,:) = G(:,:,cell_ind)';
        req_cell_ind = 1:L;
        req_cell_ind(req_cell_ind == cell_ind) = [];
        req_cell_ind(req_cell_ind == next_cell_ind) = [];
        sz = dof;
        for i = req_cell_ind
            for j = 1:K
                mat_req(sz + 1:sz + dof,:) = U(:,:,j,i)'*H(:,:,j,i,cell_ind);
                sz = sz + dof;
            end
        end
        V(:,:,user_ind,cell_ind) = null(mat_req);
    end
end
sz = size(V,2);
% R_S = 0;
R = zeros(1);
for cell_ind = 1
    H_tilda = zeros(dof,sz,K);
    for k = 1:K
        H_tilda(:,:,k) = U(:,:,k,cell_ind)'*H(:,:,k,cell_ind,cell_ind)*V(:,:,k,cell_ind);
    end
    R(cell_ind) = SumCapacityOfUserSet_DPC(1:K, H_tilda, SNR, sz, dof);
end
R_S = mean(R);


function R_S = SumCapacityOfUserSet_DPC(comb, H, SNR, n_t, n_r)

% R_S = 0;
set_card = length(comb);
P = getMACCovMat(comb, H, SNR, set_card, n_t, n_r);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Note: While solving dual MAC problem it is assumed that decoding order is
% 1,2,...,K while the eqn for BC capacity (3) in "Duality, Achievable Rates,
% and Sum-Rate Capacity of Gaussian MIMO Broadcast Channels" is written
% assuming encoding order 1,2,...,K. Therefore here after solving dual
% problem the encoding order will reverse i.e. in BC K,K-1,..,1 and the
% capacity eqn will change as below:

% Q = getMACToBCCovMat(comb, P, H, set_card, n_t, n_r);
% for j = 1:set_card
%     mat = zeros(n_t,n_t);
%     for i = 1:j-1
%         mat = mat + Q(:,:,i);
%     end
%     mat_num = mat + Q(:,:,j);
%     num = real(det(eye(n_r) + H(:,:,comb(j))*mat_num*H(:,:,comb(j))'));
%     den = real(det(eye(n_r) + H(:,:,comb(j))*mat*H(:,:,comb(j))'));
%     R_S = R_S + log2(num/den);
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Alternate way to evaluate capacity (MAC)
mat = eye(n_t);
for j = 1:set_card
    mat = mat + H(:,:,comb(j))'*P(:,:,j)*H(:,:,comb(j));
end
R_S = log2(real(det(mat)));

% function Q_bar = getMACToBCCovMat(comb, P, H, K, n_t, n_r)
% 
% Q_bar = zeros(n_t,n_t,K);
% 
% for j = 1:K
%     mat = zeros(n_t,n_t);
%     for l = 1:j-1
%         mat = mat + Q_bar(:,:,l);
%     end
%     
%     D = eye(n_r) + H(:,:,comb(j))*mat*H(:,:,comb(j))';
%     E = eye(n_t);
%     for l = j+1 : K
%         E = E + H(:,:,comb(l))' * P(:,:,l) * H(:,:,comb(l));
%     end
%     [F,~,G] = svd(E^-0.5 * H(:,:,comb(j))' * D^-0.5, 0);
%     Q_bar(:,:,j) = (E^-0.5)*F*(G')*(D^0.5)*P(:,:,j)*(D^0.5)*G*(F')*(E^-0.5);
% end


function P = getMACCovMat(comb,H, Pow, K, n_t, n_r)

G = zeros(n_r,n_t,K);
P = zeros(n_r,n_r,K);
for j = 1:K
    P(:,:,j) = eye(n_r)*Pow/(K*n_r);
end
for n = 1:75
    for j = 1:K
        mat = eye(n_t);
        for i = 1:K
            if i == j, continue;end
            mat = mat + H(:,:,comb(i))' *P(:,:,i)* H(:,:,comb(i));
        end
        G(:,:,j) = H(:,:,comb(j)) * (mat^-0.5);
    end
    S = getWaterFillCovMat(G,n_r,K, Pow);
    P_nM1 = P;
    MSE = zeros(1,K);
    for j = 1:K
        if K <= 2
            P(:,:,j) = S(:,:,j);
        else
            P(:,:,j) = 1/K * S(:,:,j) + (K-1)/K * P(:,:,j);
        end
        MSE(j) = sum(sum((P(:,:,j) - P_nM1(:,:,j)).^2))/(n_t*n_r);
    end
    if sum(MSE)/K < 1e-7, break;end
end


function S = getWaterFillCovMat(G,n_r,K, Pow)

S = zeros(n_r,n_r,K);
lambda = zeros(n_r*K,1);
for i = 1:K
    lambda((i-1)*n_r+1:i*n_r,1) = real(eig(G(:,:,i)*G(:,:,i)'));
end

pl = Water_filling(lambda, Pow);

for i = 1:K
    [U,~] = eig(G(:,:,i)*G(:,:,i)');
    Lam = diag(pl((i-1)*n_r+1:i*n_r));
    S(:,:,i) = U*Lam*U';
end


function R_F = SumCapacity_grouping_Optimal(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
K = min(floor((n_t/DoF-1)/(BS-1)), 1+floor((n_r/DoF-1)/(BS-1)));
N_sim = 1000;
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            comb = combntns(1:nusers,K);
            all_perm = allcomb_gg(1:size(comb,1), BS);
            for k = 1:size(all_perm,1)
                user_comb = comb(all_perm(k,:)',:);
                R_S = Sum_capacity_grouping(H, SNR_vect(SNR_ind), user_comb, K,BS, DoF);
                if R_S > R_final(avg_sim_ind)
                    R_final(avg_sim_ind) = R_S;
                end
            end
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = Sum_capacity_MIMO_IC(H, SNR,comb)

global n_t
global n_r

U = zeros(n_r,n_r,2);
V = zeros(n_t,n_t,2);

[U(:,:,1),~,V(:,:,1)] = svd(H(:,:,2,1,comb(1)));
[U(:,:,2),~,V(:,:,2)] = svd(H(:,:,1,2,comb(2)));
% For H_1 we will take the last N columns
H_1_tmp = U(:,:,1)'*H(:,:,1,1,comb(1))*V(:,:,2);
H_1 = H_1_tmp(:,n_t-n_r+1:end);
% For H_2 we will take the first min(M-N,N) inputs
H_2_tmp = U(:,:,2)'*H(:,:,2,2,comb(2))*V(:,:,1);
H_2 = H_2_tmp(1:min(n_t-n_r,n_r),n_r+1:end);

Pow = 10^(SNR / 10);
lambda1 = eig(H_1*H_1');
pl1 = Water_filling(lambda1, Pow);
lambda2 = eig(H_2*H_2');
pl2 = Water_filling(lambda2, Pow);
R_F = sum(log2(1 + lambda1.*pl1)) + sum(log2(1 + lambda2.*pl2));


function R_F = SumCapacity_grouping(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        comb = repmat(1:nusers,BS,1);
        for avg_sim_ind = 1:N_sim

%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>         
            H = (randn(n_r,n_t,nusers,BS,BS) + 1i*randn(n_r,n_t,nusers,BS,BS))/sqrt(2);
            R_final(avg_sim_ind) = Sum_capacity_grouping(H, SNR_vect(SNR_ind),comb, nusers, BS, DoF);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_three_cell_alignment(SNR_vect, BS, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        comb = repmat(1:nusers,BS,1);
        for avg_sim_ind = 1:N_sim

%             H = (randn(n_r,n_t,3,3,nusers) + 1i*randn(n_r,n_t,3,3,nusers))/sqrt(2);
            H = (randn(n_r,n_t,3,3,nusers) + 1i*randn(n_r,n_t,3,3,nusers))/sqrt(2);
%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>         
            R_final(avg_sim_ind) = Sum_capacity_three_cell(H, SNR_vect(SNR_ind),comb, nusers, DoF);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function [R_F,input_approach] = SumCapacity_cyclic_model2(SNR_vect, BS, nusers)

global n_t
global n_r
R_F = zeros(length(SNR_vect), 1);
time = zeros(length(SNR_vect), 1);
input_approach = [];

fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end

fprintf('Enter the type of approach in Model 2, curretly available [b, c] :\n');
input_approach = input('','s');
err_flag = 0;
switch input_approach
    
    case 'b'
        if n_t < nusers*n_r || n_r < (nusers+1)*DoF
            err_flag = 1;
        end
        Sum_capacity = @ Sum_capacity_cyclic_model2_b;
        
    case 'c'
        if n_t < 2*nusers*DoF || n_r < (nusers+1)*DoF
            err_flag = 1;
        end
        Sum_capacity = @ Sum_capacity_cyclic_model2_c;
        
    otherwise
        fprintf('Input approach doesn''t exist, program will terminate!\n');
        return
end
if err_flag
    fprintf('Antenna values doesn''t met constraints, program will terminate!\n');
    return
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
%         eval(['load ' str]);
    dlg.FractionComplete = 0;
    dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
    tic
    comb = repmat(1:nusers,BS,1);
    for avg_sim_ind = 1:N_sim

        H = (randn(n_r,n_t,BS,BS,nusers) + 1i*randn(n_r,n_t,BS,BS,nusers))/sqrt(2);
%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>         
        R_final(avg_sim_ind) = Sum_capacity(H, SNR_vect(SNR_ind),comb, BS, nusers, DoF);
        if ~rem(avg_sim_ind,10)
            dlg.FractionComplete = avg_sim_ind / N_sim;
            dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
        end
    end 
    R_F(SNR_ind) = mean(R_final);
    time(SNR_ind) = toc;
    fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(time));
end


function [R_F,input_approach] = SumCapacity_cyclic_model2_Capacity(SNR_vect, L, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
input_approach = [];

fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end

fprintf('Enter the type of approach in Model 2, curretly available [b, c] :\n');
input_approach = input('','s');
% err_flag = 0;
switch input_approach
    
    case 'b'
%         if n_t < users_num*n_r || n_r < (users_num+1)*DoF,
%             err_flag = 1;
%         end
        K = n_t / n_r;
        Sum_capacity = @ Sum_capacity_cyclic_model2_b;
        
    case 'c'
%         if n_t < 2*users_num*DoF || n_r < (users_num+1)*DoF,
%             err_flag = 1;
%         end
        K = n_t / 2 / DoF;
        Sum_capacity = @ Sum_capacity_cyclic_model2_c;
        
    otherwise
        fprintf('Input approach doesn''t exist, program will terminate!\n');
        return
end
% if err_flag,
%     fprintf('Antenna values doesn''t met constraints, program will terminate!\n');
%     return
% end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        comb = repmat(1:nusers,L,1);
        for avg_sim_ind = 1:N_sim

            H = (randn(n_r,n_t,L,L,nusers) + 1i*randn(n_r,n_t,L,L,nusers))/sqrt(2);
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('new',H(:,:,cell_ind,cell_ind,:), K, 1:nusers);
            end
           
            R_temp = 0;
            maxValInd_exist = 0;
            for cell_ind = 1:L
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
%                     T_set(T_set == other_users) = [];
                    T_set(other_users) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        R_S = Sum_capacity(H, SNR_vect(SNR_ind),comb, L, K, DoF);
                        if R_S > R_temp
                            R_temp = R_S;
                            maxValInd = user_ind;
                            maxValInd_exist = 1;
                        end
                    end
                    if maxValInd_exist
                        user_comb(cell_ind,k) = maxValInd;
                    end
                    maxValInd_exist = 0;
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_three_cell_Chordal(SNR_vect, L, users_num)

global n_t
global n_r

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
% DoF = 2;
% K = min(floor((n_t/DoF-1)/(L-1)), 1+floor((n_r/DoF-1)/(L-1)));
K = 2;
if K<=1
    fprintf('Unsupported value input (K<=1), program will terminate!\n');
    return
end 
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = (randn(n_r,n_t,L,L,nusers) + 1i*randn(n_r,n_t,L,L,nusers))/sqrt(2);
%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('new',H(:,:,cell_ind,cell_ind,:), K, 1:nusers);
            end

            R_temp = 0;
            for cell_ind = 1:L
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
%                     T_set(T_set == other_users) = [];
                    T_set(other_users) = [];
                    comb_temp = user_comb;
                    Ent_temp = 0;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        Ent = compute_three_cell_chordal_dist(H, comb_temp, K, DoF, [k,cell_ind]);
                        if Ent > Ent_temp
                            Ent_temp = Ent;
                            maxValInd = user_ind;
                        end
                    end
                    comb_temp = user_comb;
                    comb_temp(cell_ind,k) = maxValInd;
                    R_S = Sum_capacity_three_cell(H, SNR_vect(SNR_ind),comb_temp, K, DoF);
                    if R_S > R_temp
                        user_comb = comb_temp;
                        R_temp = R_S;
                    end
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_grouping_Chordal(SNR_vect, L, users_num)

global n_t
global n_r

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
% DoF = 2;
K = min(floor((n_t/DoF-1)/(L-1)), 1+floor((n_r/DoF-1)/(L-1)));
K = 2;
if K<=1
    fprintf('Unsupported value input (K<=1), program will terminate!\n');
    return
end 
N_sim = 1000;
% alpha = 0:0.1:1;
% alp_ind = 1;
% R_F = zeros(1,length(alpha));
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
%         for alp_ind = 1:length(alpha)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
        eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        for avg_sim_ind = 1:N_sim

%             H = (randn(n_r,n_t,nusers,L,L) + 1i*randn(n_r,n_t,nusers,L,L))/sqrt(2);
            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('old',H(:,:,:,cell_ind,cell_ind), K, 1:nusers);
            end
            [U,G] = getInitialReceiverMatrices(H, user_comb, K, L, DoF);

            R_temp = 0;
            for cell_ind = 1:L
                U_change_flag = 0;
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
                    T_set(T_set == other_users) = [];
                    comb_temp = user_comb;
                    Ent_temp = 0;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        [U_current,G_current] = UpdateUandG_Matrices(H,comb_temp, U, G, K, L, DoF, [k,cell_ind]);
                        Ent = compute_chordal_dist(H, comb_temp, K, L, DoF, [k,cell_ind], U_current, G_current);
                        if Ent > Ent_temp
                            Ent_temp = Ent;
                            maxValInd = user_ind;
                            maxValuser_U = U_current;
                            maxValuser_G = G_current;
                        end
                    end
                    comb_temp = user_comb;
                    comb_temp(cell_ind,k) = maxValInd;
                    R_S = Sum_capacity_grouping_mod(H, SNR_vect(SNR_ind), comb_temp, K,L, DoF, maxValuser_U,maxValuser_G);
                    if R_S > R_temp
                        user_comb = comb_temp;
                        R_temp = R_S;
                        U_req = maxValuser_U;
                        G_req = maxValuser_G;
                        U_change_flag = 1;
                    end
                end
                if U_change_flag
                    U = U_req;
                    G = G_req;
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function [U,G] = getInitialReceiverMatrices(H, comb, K, L, dof)

global n_t
global n_r

U = zeros(n_r,dof,K,L);
G = zeros(n_t,dof,L);
% W = zeros(dof,dof,K,L);

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if cell_ind + 1 > L, next_cell_ind = 1;end
    if K <= 3
        mat_req = zeros(K*n_t,n_t + K*n_r);
        sz_eye = 1;
        sz_H_colm = n_t;
        sz_H_row = 1;
        for i = 1:K
            mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
            mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,comb(next_cell_ind,i),next_cell_ind,cell_ind)';
            sz_eye = sz_eye + n_t;
            sz_H_row = sz_H_row + n_t;
            sz_H_colm = sz_H_colm + n_r;
        end
        X = null(mat_req);
        G(:,:,cell_ind) = X(1:n_t,1:dof);
        sz = n_t;
        for i = 1:K
            U(:,:,i,next_cell_ind) = X(sz+1:sz+n_r,1:dof);
            sz = sz + n_r;
        end
    else
        G_tilda = zeros(n_t,n_r,K);
        U_tilda = zeros(n_r,n_r,K);
        for k = 1:K
            mat = [eye(n_t) -H(:,:,comb(next_cell_ind,k),next_cell_ind,cell_ind)'];
            X_tilda = null(mat);
            G_tilda(:,:,k) = X_tilda(1:n_t,:);
            U_tilda(:,:,k) = X_tilda(n_t+1:end,:);
        end
        [G_get,U_get] = Compute_common_subspace(G_tilda, U_tilda);
        G(:,:,cell_ind) = G_get;
        U(:,:,:,next_cell_ind) = U_get;
    end
%     for i = 1:K
%         W(:,:,i,next_cell_ind) = inv(sqrtm(U(:,:,i,next_cell_ind)'*U(:,:,i,next_cell_ind)));
%     end
end


% function [U_current,G_current,W_current] = UpdateUandG_Matrices(H, comb, U, G, W, K, L, dof, user_id)
function [U_current,G_current] = UpdateUandG_Matrices(H, comb, U, G, K, L, dof, user_id)

global n_t
global n_r

U_current = U;
G_current = G;
% W_current = W;
prev_cell_ind = user_id(2) - 1;
if user_id(2) == 1, prev_cell_ind = L;end

if K <= 3
    mat_req = zeros(K*n_t,n_t + K*n_r);
    sz_eye = 1;
    sz_H_colm = n_t;
    sz_H_row = 1;
    for i = 1:K
        mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
        mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,comb(user_id(2),i),user_id(2),prev_cell_ind)';
        sz_eye = sz_eye + n_t;
        sz_H_row = sz_H_row + n_t;
        sz_H_colm = sz_H_colm + n_r;
    end
    X = null(mat_req);
    G_current(:,:,prev_cell_ind) = X(1:n_t,1:dof);
    sz = n_t;
    for i = 1:K
        U_current(:,:,i,user_id(2)) = X(sz+1:sz+n_r,1:dof);
        sz = sz + n_r;
    end
else
    G_tilda = zeros(n_t,n_r,K);
    U_tilda = zeros(n_r,n_r,K);
    for k = 1:K
        mat = [eye(n_t) -H(:,:,comb(user_id(2),k),user_id(2),prev_cell_ind)'];
        X_tilda = null(mat);
        G_tilda(:,:,k) = X_tilda(1:n_t,:);
        U_tilda(:,:,k) = X_tilda(n_t+1:end,:);
    end
    [G_get,U_get] = Compute_common_subspace(G_tilda, U_tilda);
    U_current(:,:,:,user_id(2)) = U_get;
    G_current(:,:,prev_cell_ind) = G_get;
end
% for i = 1:K
%     W_current(:,:,i,user_id(2)) = inv(sqrtm(U_current(:,:,i,user_id(2))'*U_current(:,:,i,user_id(2))));
% end


function met_val = compute_three_cell_chordal_dist(H, comb, K, dof, user_id)

global n_t
global n_r

G = zeros(n_t, n_t, 3, 3);
W_bar = zeros(n_t, dof, 3);
W = zeros(n_r, dof, 3, K);

for cell_ind = 1:3
    for other_cell_ind = 1:3
        if other_cell_ind == cell_ind, continue;end
        sz = 1;
        for user_ind = 1:K
            G(:,sz:sz+n_r-1,cell_ind,other_cell_ind) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,user_ind))';
            sz = sz + n_r;
        end
    end
end
E = (G(:,:,3,1) \ G(:,:,3,2)) * (G(:,:,1,2) \ G(:,:,1,3)) * (G(:,:,2,3) \ G(:,:,2,1));
F = G(:,:,3,2) \ G(:,:,3,1);
C = G(:,:,2,3) \ G(:,:,2,1);
[Eig_vect,~] = eig(E);
W_bar(:,:,1) = Eig_vect(:,1:dof);
W_bar(:,:,2) = F * W_bar(:,:,1);
W_bar(:,:,3) = C * W_bar(:,:,1);

for cell_ind = 1:3
    sz = 1;
    for user_ind = 1:K
        W(:,:,cell_ind,user_ind) = W_bar(sz:sz+n_r-1,:,cell_ind);
        sz = sz + n_r;
    end
end

H_A = zeros(n_t, K*dof);
sz = 1;
for k = 1:K
    H_A(:,sz:sz + dof-1) = H(:,:,user_id(2),user_id(2),comb(user_id(2),k))' * W(:,:,user_id(2),k);
    sz = sz + dof;
end
% H_A = H(:,:,user_id(2),user_id(2),comb(user_id(2),user_id(1)))' * W(:,:,user_id(2),user_id(1));
% H_B_hor = zeros(n_t,(3*K-1)*dof);
H_B_hor = zeros(n_t,(2*K)*dof);
sz_hor = 1;
for cell_ind = 1:3
    if cell_ind == user_id(2),continue;end
    for user_ind = 1:K
%         if user_ind == user_id(1) && cell_ind == user_id(2),continue;end
        H_B_hor(:,sz_hor:sz_hor+dof-1) = H(:,:,user_id(2),cell_ind,comb(cell_ind,user_ind))'*...
            W(:,:,cell_ind,user_ind);
        sz_hor = sz_hor + dof;
    end
end

A = GSO(H_A);
B = GSO(H_B_hor);
B_req = B(:,1:(2*K-1)*dof);
met_val = norm((A*A' - B_req*B_req'),'fro');


function [met_val] = compute_chordal_dist(H, comb, K, L, dof, user_id, U_current, G_current)

global n_t
% global n_r

next_cell_ind = user_id(2) + 1;
if user_id(2) + 1 > L, next_cell_ind = 1;end

H_A = H(:,:,comb(user_id(2),user_id(1)),user_id(2),user_id(2))' * U_current(:,:,user_id(1),user_id(2));
% H_A = H(:,:,comb(user_id(2),user_id(1)),user_id(2),user_id(2))';
cell_set = 1:L;
cell_set(cell_set == next_cell_ind) = [];
H_B_hor = zeros(n_t,(K*(L-1)-1)*dof);
% H_B_hor = zeros(n_t,(K*(L-1)-1)*n_r);
sz_hor = 1;
for cell_ind = cell_set
    for user_ind = 1:K
        if user_ind == user_id(1) && cell_ind == user_id(2),continue;end
        H_B_hor(:,sz_hor:sz_hor+dof-1) = H(:,:,comb(cell_ind,user_ind),cell_ind,user_id(2))'*...
            U_current(:,:,user_ind,cell_ind);
        sz_hor = sz_hor + dof;
%         H_B_hor(:,sz_hor:sz_hor+n_r-1) = H(:,:,comb(cell_ind,user_ind),cell_ind,user_id(2))';
%         sz_hor = sz_hor + n_r;
    end
end

% H_B_hor = [H_B_hor,H(:,:,comb(next_cell_ind,1),next_cell_ind,user_id(2))'*...
%             U_current(:,:,1,next_cell_ind)];
H_B_hor = [H_B_hor,G_current(:,:,user_id(2))];
% H_B_hor = [H_B_hor,H(:,:,comb(next_cell_ind,1),next_cell_ind,user_id(2))'];

met_val = cho_dist(H_A,H_B_hor);


function R_F = SumCapacity_grouping_Capacity(SNR_vect, L, users_num)

global n_t
global n_r
% global data_gg

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
% DoF = 2;
K = min(floor((n_t/DoF-1)/(L-1)), 1+floor((n_r/DoF-1)/(L-1)));
if K<=1
    fprintf('Unsupported value input (K<=1), program will terminate!\n');
    return
end 
N_sim = 1000;
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
        eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        for avg_sim_ind = 1:N_sim
            
%             H = (randn(n_r,n_t,nusers,L,L) + sqrt(-1)*randn(n_r,n_t,nusers,L,L))/sqrt(2);

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('old',H(:,:,:,cell_ind,cell_ind), K, 1:nusers);
            end
            
            [U,G] = getInitialReceiverMatrices(H, user_comb, K, L, DoF);
            R_temp = 0;
            maxValInd_exist = 0;
            for cell_ind = 1:L
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
%                     T_set(T_set == other_users) = [];
                    T_set(other_users) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        [U_current,G_current] = UpdateUandG_Matrices(H,comb_temp, U, G, K, L, DoF, [k,cell_ind]);
                        R_S = Sum_capacity_grouping_mod(H, SNR_vect(SNR_ind), comb_temp, K, L, DoF, U_current, G_current);

                        if R_S > R_temp
                            R_temp = R_S;
                            maxValInd = user_ind;
                            maxValInd_exist = 1;
                            maxValUser_UandG = {U_current,G_current};
                        end
                    end
                    if maxValInd_exist
                        user_comb(cell_ind,k) = maxValInd;
                        U = maxValUser_UandG{1};
                        G = maxValUser_UandG{2};
                        UandG_flag = 1;
                    end
                    maxValInd_exist = 0;
                end
                if UandG_flag
                    U = maxValUser_UandG{1};
                    G = maxValUser_UandG{2};
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_three_cell_Capacity(SNR_vect, L, users_num)

global n_t
global n_r

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end

K = 2;
N_sim = 1000;
% dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
%         eval(['load ' str]);
%         dlg.FractionComplete = 0;
%         dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        progress(0);
        tic
        for avg_sim_ind = 1:N_sim
            
%             H = (randn(n_r,n_t,nusers,L,L) + sqrt(-1)*randn(n_r,n_t,nusers,L,L))/sqrt(2);
            H = (randn(n_r,n_t,L,L,nusers) + sqrt(-1)*randn(n_r,n_t,L,L,nusers))/sqrt(2);

%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('new',H(:,:,cell_ind,cell_ind,:), K, 1:nusers);
            end
            
            R_temp = 0;
            maxValInd_exist = 0;
            for cell_ind = 1:L
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
%                     T_set(T_set == other_users) = [];
                    T_set(other_users) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        R_S = Sum_capacity_three_cell(H, SNR_vect(SNR_ind), comb_temp, K, DoF);

                        if R_S > R_temp
                            R_temp = R_S;
                            maxValInd = user_ind;
                            maxValInd_exist = 1;
                        end
                    end
                    if maxValInd_exist
                        user_comb(cell_ind,k) = maxValInd;
                    end
                    maxValInd_exist = 0;
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
%                 dlg.FractionComplete = avg_sim_ind / N_sim;
%                 dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
                progress(avg_sim_ind / N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_grouping_CapacityTest(SNR_vect, L, users_num)

global n_t
global n_r

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
% fprintf('Enter the value of DoF for each user in each cell :\n');
% DoF_str = input('','s');
% [DoF,status] = str2num(DoF_str);
% if ~status
%     fprintf('Wrong value input, program will terminate!\n');
%     return
% end
DoF = 2;
K = min(floor((n_t/DoF-1)/(L-1)), 1+floor((n_r/DoF-1)/(L-1)));
if K<1
    fprintf('Unsupported value input (K<=1), program will terminate!\n');
    return
end 
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            
            user_comb = zeros(L,2);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('old',H(:,:,:,cell_ind,cell_ind), 2, 1:nusers);
            end
            
            [U,G] = getInitialReceiverMatrices(H, user_comb, 2, L, DoF);
            R_temp = 0;
            maxValInd_exist = 0;
            for cell_ind = 1:L
                user_set = 1:nusers;
                for k = 1:2
                    T_set = user_set;
                    other_user_ind = 1:2;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
                    T_set(T_set == other_users) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        [U_current,G_current] = UpdateUandG_Matrices(H,comb_temp, U, G, 2, L, DoF, [k,cell_ind]);
                        R_S = Sum_capacity_grouping_mod(H, SNR_vect(SNR_ind), comb_temp, 2, L, DoF, U_current, G_current);

                        if R_S > R_temp
                            R_temp = R_S;
                            maxValInd = user_ind;
                            maxValInd_exist = 1;
                            maxValUser_UandG = {U_current,G_current};
                        end
                    end
                    if maxValInd_exist,
                        user_comb(cell_ind,k) = maxValInd;
                        U = maxValUser_UandG{1};
                        G = maxValUser_UandG{2};
                    end
                    maxValInd_exist = 0;
                end
            end
            
            for k = 3:K
                for l = 1:L
                    T_set = 1:nusers;
                    T_set(user_comb(l,1:k-1)) = [];
                    user_comb(l,k) = getInitialUserSet('old',H(:,:,:,l,l), 1, T_set);
                end
                user_set = 1:nusers;
                for cell_ind = 1:L
                    T_set = user_set;
%                     other_user_ind = 1:K-1;
%                     other_user_ind(other_user_ind == k) = [];
%                     other_users = user_comb(cell_ind,other_user_ind);
                    T_set(user_comb(cell_ind,1:K-1)) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,end) = user_ind;
%                         [U_current,G_current] = UpdateUandG_Matrices(H,comb_temp, U, G, k, L, DoF, [k,cell_ind]);
                        [U,G] = getInitialReceiverMatrices(H, comb_temp, k, L, DoF);
                        R_S = Sum_capacity_grouping_mod(H, SNR_vect(SNR_ind), comb_temp, k, L, DoF, U, G);

                        if R_S > R_temp
                            R_temp = R_S;
                            maxValInd = user_ind;
                            maxValInd_exist = 1;
%                             maxValUser_UandG = {U_current,G_current};
                        end
                    end
                    if maxValInd_exist
                        user_comb(cell_ind,end) = maxValInd;
%                         U = maxValUser_UandG{1};
%                         G = maxValUser_UandG{2};
                    end
                    maxValInd_exist = 0;
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_S = Sum_capacity_grouping_exp(H, SNR, comb, K, L, dof, user_id)

global n_t
global n_r
global data_gg

SNR = 10^(SNR/10);
U = zeros(n_r,dof,K,L);
V = zeros(n_t,dof,K,L);
G = zeros(n_t,dof,L);

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if cell_ind + 1 > L, next_cell_ind = 1;end
    if K <= 3
        mat_req = zeros(K*n_t,n_t + K*n_r);
        sz_eye = 1;
        sz_H_colm = n_t;
        sz_H_row = 1;
        for i = 1:K
            mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
            mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,comb(next_cell_ind,i),next_cell_ind,cell_ind)';
            sz_eye = sz_eye + n_t;
            sz_H_row = sz_H_row + n_t;
            sz_H_colm = sz_H_colm + n_r;
        end
        X = null(mat_req);
        G(:,:,cell_ind) = X(1:n_t,1:dof);
        sz = n_t;
        for i = 1:K
            U(:,:,i,next_cell_ind) = X(sz+1:sz+n_r,1:dof);
            sz = sz + n_r;
        end
    else
        G_tilda = zeros(n_t,n_r,K);
        U_tilda = zeros(n_r,n_r,K);
        for k = 1:K
            mat = [eye(n_t) -H(:,:,comb(next_cell_ind,k),next_cell_ind,cell_ind)'];
            X_tilda = null(mat);
            G_tilda(:,:,k) = X_tilda(1:n_t,:);
            U_tilda(:,:,k) = X_tilda(n_t+1:end,:);
        end
        [G_get,U_get] = Compute_common_subspace(G_tilda, U_tilda);
        G(:,:,cell_ind) = G_get;
        U(:,:,:,next_cell_ind) = U_get;
    end
end

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if next_cell_ind > L, next_cell_ind = 1;end
    for user_ind = 1:K
        mat_req = zeros(K*(L-1)*dof,n_t);
        mat_req(1:dof,:) = G(:,:,cell_ind)';
        req_cell_ind = 1:L;
        req_cell_ind(req_cell_ind == cell_ind) = [];
        req_cell_ind(req_cell_ind == next_cell_ind) = [];
        sz = dof;
        for i = req_cell_ind
            for j = 1:K
                mat_req(sz + 1:sz + dof,:) = U(:,:,j,i)'*H(:,:,comb(i,j),i,cell_ind);
                sz = sz + dof;
            end
        end
        req_user_ind = 1:K;
        req_user_ind(req_user_ind == user_ind) = [];
        for j = req_user_ind
            mat_req(sz + 1:sz + dof,:) = U(:,:,j,cell_ind)'*H(:,:,comb(cell_ind,j),cell_ind,cell_ind);
            sz = sz + dof;
        end
        V(:,:,user_ind,cell_ind) = null(mat_req);
    end
end

R_S = 0;
Pow = SNR/dof/K; % Assume equal power allocation and V is orthonormal
for cell_ind = 1:L
    for user_ind = 1:K
%         R_S = R_S + Compute_SINR_Capacity(H, U, V, SNR, user_ind, cell_ind, K, L, dof);
        mat = U(:,:,user_ind,cell_ind)'*H(:,:,comb(cell_ind,user_ind),cell_ind,cell_ind)*V(:,:,user_ind,cell_ind);
        R_S = R_S + log2(real(det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind) + Pow*(mat*mat'))/...
            det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind))));
    end
end

H_A = H(:,:,comb(user_id(2),user_id(1)),user_id(2),user_id(2))' * U(:,:,user_id(1),user_id(2));
cell_set = 1:L;
cell_set(cell_set == next_cell_ind) = [];
H_B_hor = zeros(n_t,(K*(L-1)-1)*dof);
sz_hor = 1;
for cell_ind = cell_set
    for user_ind = 1:K
        if user_ind == user_id(1) && cell_ind == user_id(2),continue;end
        H_B_hor(:,sz_hor:sz_hor+dof-1) = H(:,:,comb(cell_ind,user_ind),cell_ind,user_id(2))'*...
            U(:,:,user_ind,cell_ind);
        sz_hor = sz_hor + dof;
    end
end

H_B_hor = [H_B_hor,H(:,:,comb(next_cell_ind,1),next_cell_ind,user_id(2))'*...
            U(:,:,1,next_cell_ind)];
% alpha = -2;
met_val = cho_dist(H_A,H_B_hor);
H_req = V(:,:,user_id(1),user_id(2))'*H_A;

data_gg = [data_gg ; {R_S}, {met_val}, {norm(H_A)}, {norm(H_req)}];


function user_comb = getInitialUserSet(dir, H, K, tot_users)


Frob_norm = zeros(1,length(tot_users));
if strcmp(dir, 'new')
    for user_ind = tot_users
        Frob_norm(user_ind) = norm(H(:,:,:,:,user_ind),'fro');
    end
    
elseif strcmp(dir, 'old')
    for user_ind = tot_users
        Frob_norm(user_ind) = norm(H(:,:,user_ind),'fro');
    end
end

[~,sorted_ind] = sort(Frob_norm,'descend');
user_comb = sorted_ind(1:K);

% FrobTemp = 0;
% O_set = 1:nusers;
% for i = 1:nusers
%     FrobVal = trace(H(:,:,i) * H(:,:,i)');
%     if FrobVal > FrobTemp,
%         FrobTemp = FrobVal;
%         maxValInd = i;
%     end
% end
% O_set(O_set == maxValInd) = [];
% G_set = maxValInd;
% Q = GSO(H(:,:,G_set(1)).');
% V = Q.';
% 
% for i = 2:K
%     F_Temp = 0;
%     for k = O_set
%         Hk_tilda = H(:,:,k) - H(:,:,k) * (V') * V;
%         FrobNorm = norm(Hk_tilda, 'fro')^2;
%         H_cap = zeros(n_r*(i-1), n_t);
%         for s = G_set
%             count = 1;
%             for g = G_set
%                 if g == s, continue;end
%                 H_cap(n_r*(count-1)+1:n_r*count,:) = H(:,:,g);
%                 count = count + 1;
%             end
%             H_cap(end-n_r+1:end,:) = H(:,:,k);
%             Q = GSO(H_cap.');
%             W = Q.';
%             Hs_tilda = H(:,:,s) - H(:,:,s) * (W') * W;
%             FrobNorm = FrobNorm + norm(Hs_tilda, 'fro')^2;
%         end
%         if FrobNorm > F_Temp
%             F_Temp = FrobNorm;
%             maxValInd = k;
%         end
%     end
%     O_set(O_set == maxValInd) = [];
%     G_set = [G_set,maxValInd];
%     Q = GSO(H(:,:,maxValInd).');
%     Vs = Q.';
%     V = [V;Vs];
% end
% user_comb = G_set;

% Pow = 10 ^ (SNR/10);
% Ent = zeros(1,nusers);
% for user_ind = 1:nusers
%     Ent(user_ind) = real(det(eye(n_r) + Pow*H(:,:,user_ind)*...
%         H(:,:,user_ind)'));
% end
% [~,sorted_ind] = sort(Ent,'descend');
% user_comb = sorted_ind(1:K);
% 
% Ent_temp = 0;
% user_set = 1:nusers;
% for user_ind = user_set
%     Ent = det(eye(n_r) + Pow/n_t*H(:,:,user_ind)*H(:,:,user_ind)');
%     if Ent > Ent_temp,
%         Ent_temp = Ent;
%         maxValInd = user_ind;
%     end
% end
% CEnt_temp = 0;
% u1_ind = maxValInd;
% user_set(user_set == maxValInd) = [];
% for user_ind = user_set
%     H12 = eye(n_t) + Pow/n_t*H(:,:,u1_ind)'*H(:,:,u1_ind) + Pow*H(:,:,user_ind)'*...
%         H(:,:,user_ind);
%     H2 = eye(n_t) + Pow/n_t*H(:,:,user_ind)'*H(:,:,user_ind);
%     CEnt = (det(H12))^2 / det(H2);
%     if CEnt > CEnt_temp,
%         CEnt_temp = CEnt;
%         maxValInd = user_ind;
%     end
% end
% user_comb = [u1_ind,maxValInd];


function R_S = Sum_capacity_grouping_mod(H, SNR, comb, K, L, dof, U, G)

global n_t
% global n_r

SNR = 10^(SNR/10);
V = zeros(n_t,dof,K,L);


for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if next_cell_ind > L, next_cell_ind = 1;end
    for user_ind = 1:K
        mat_req = zeros(K*(L-1)*dof,n_t);
        mat_req(1:dof,:) = G(:,:,cell_ind)';
        req_cell_ind = 1:L;
        req_cell_ind(req_cell_ind == cell_ind) = [];
        req_cell_ind(req_cell_ind == next_cell_ind) = [];
        sz = dof;
        for i = req_cell_ind
            for j = 1:K
                mat_req(sz + 1:sz + dof,:) = U(:,:,j,i)'*H(:,:,comb(i,j),i,cell_ind);
                sz = sz + dof;
            end
        end
        req_user_ind = 1:K;
        req_user_ind(req_user_ind == user_ind) = [];
        for j = req_user_ind
            mat_req(sz + 1:sz + dof,:) = U(:,:,j,cell_ind)'*H(:,:,comb(cell_ind,j),cell_ind,cell_ind);
            sz = sz + dof;
        end
        X = null(mat_req);
        V(:,:,user_ind,cell_ind) = X(:,1:dof);
    end
end


R_S = 0;
Pow = SNR/dof/K; % Assume equal power allocation and V is orthonormal
for cell_ind = 1:L
    for user_ind = 1:K
%         R_S = R_S + Compute_SINR_Capacity(H, U, V, SNR, user_ind, cell_ind, K, L, dof);
        mat = U(:,:,user_ind,cell_ind)'*H(:,:,comb(cell_ind,user_ind),cell_ind,cell_ind)*V(:,:,user_ind,cell_ind);
        R_S = R_S + log2(real(det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind) + Pow*(mat*mat'))/...
            det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind))));
    end
end

% R_S = 0;
% for cell_ind = 1:L
%     lambda = zeros(dof*K,1);
%     sz = 1;
%     for user_ind = 1:K
%         W(:,:,user_ind,cell_ind) = inv(sqrtm(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind)));
%         H_eq = W(:,:,user_ind,cell_ind)*U(:,:,user_ind,cell_ind)'*H(:,:,comb(cell_ind,user_ind),cell_ind,cell_ind)*V(:,:,user_ind,cell_ind);
%         lambda(sz:sz+dof-1,:) = eig(H_eq*H_eq');
%         sz = sz + dof;
%     end
%     pl = Water_filling(lambda, SNR);
%     R_S = R_S + sum(log2(1 + pl.*lambda));
% end


function R_S = Sum_capacity_grouping(H, SNR, comb, K, L, dof)

global n_t
global n_r

SNR = 10^(SNR/10);
U = zeros(n_r,dof,K,L);
V = zeros(n_t,dof,K,L);
G = zeros(n_t,dof,L);

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if cell_ind + 1 > L, next_cell_ind = 1;end
    if K <= 3
        mat_req = zeros(K*n_t,n_t + K*n_r);
        sz_eye = 1;
        sz_H_colm = n_t;
        sz_H_row = 1;
        for i = 1:K
            mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
            mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,comb(next_cell_ind,i),next_cell_ind,cell_ind)';
            sz_eye = sz_eye + n_t;
            sz_H_row = sz_H_row + n_t;
            sz_H_colm = sz_H_colm + n_r;
        end
        X = null(mat_req);
        G(:,:,cell_ind) = X(1:n_t,1:dof);
        sz = n_t;
        for i = 1:K
            U(:,:,i,next_cell_ind) = X(sz+1:sz+n_r,1:dof);
            sz = sz + n_r;
        end
    else
        G_tilda = zeros(n_t,n_r,K);
        U_tilda = zeros(n_r,n_r,K);
        for k = 1:K
            mat = [eye(n_t) -H(:,:,comb(next_cell_ind,k),next_cell_ind,cell_ind)'];
            X_tilda = null(mat);
            G_tilda(:,:,k) = X_tilda(1:n_t,:);
            U_tilda(:,:,k) = X_tilda(n_t+1:end,:);
        end
        [G_get,U_get] = Compute_common_subspace(G_tilda, U_tilda);
        G(:,:,cell_ind) = G_get;
        U(:,:,:,next_cell_ind) = U_get;
    end
end

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if next_cell_ind > L, next_cell_ind = 1;end
    for user_ind = 1:K
        mat_req = zeros(K*(L-1)*dof,n_t);
        mat_req(1:dof,:) = G(:,:,cell_ind)';
        req_cell_ind = 1:L;
        req_cell_ind(req_cell_ind == cell_ind) = [];
        req_cell_ind(req_cell_ind == next_cell_ind) = [];
        sz = dof;
        for i = req_cell_ind
            for j = 1:K
                mat_req(sz + 1:sz + dof,:) = U(:,:,j,i)'*H(:,:,comb(i,j),i,cell_ind);
                sz = sz + dof;
            end
        end
        req_user_ind = 1:K;
        req_user_ind(req_user_ind == user_ind) = [];
        for j = req_user_ind
            mat_req(sz + 1:sz + dof,:) = U(:,:,j,cell_ind)'*H(:,:,comb(cell_ind,j),cell_ind,cell_ind);
            sz = sz + dof;
        end
        null_mat = null(mat_req);
        V(:,:,user_ind,cell_ind) = null_mat(:,1:dof);
    end
end

R_S = 0;
Pow = SNR/dof/K; % Assume equal power allocation and V is orthonormal
for cell_ind = 1:L
    for user_ind = 1:K
%         R_S = R_S + Compute_SINR_Capacity(H, U, V, SNR, user_ind, cell_ind, K, L, dof);
        mat = U(:,:,user_ind,cell_ind)'*H(:,:,comb(cell_ind,user_ind),cell_ind,cell_ind)*V(:,:,user_ind,cell_ind);
        R_S = R_S + log2(real(det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind) + Pow*(mat*mat'))/...
            det(U(:,:,user_ind,cell_ind)'*U(:,:,user_ind,cell_ind))));
    end
end


function [G,U] = Compute_common_subspace(G_tilda, U_tilda)

global n_r

sz_G = size(G_tilda);
K = size(G_tilda,3);
if K == 1
    G = G_tilda;
    U = U_tilda;
    return;
end
sz_U = size(U_tilda);
U = zeros(1,1,sz_U(3));

if ~rem(K,2)
    G = zeros(1,1,K/2);
    for i = 1:K/2
        mat = [G_tilda(:,:,i) -G_tilda(:,:,K/2+i)];
        X_tilda = null(mat);
        C1 = X_tilda(1:sz_G(2),:);
        C2 = X_tilda(sz_G(2)+1:end,:);
        null_dim = size(X_tilda, 2);
        G(1:sz_G(1),1:null_dim,i) = G_tilda(:,:,i)*C1;
        ii = i:K:sz_U(3);
        jj =  i+K/2:K:sz_U(3);
        for ind = ii
            U(1:n_r,1:null_dim,ind) = U_tilda(:,:,ind)*C1;
        end
        for ind = jj
            U(1:n_r,1:null_dim,ind) = U_tilda(:,:,ind)*C2;
        end
    end
else
    Kby2 = ceil(K/2);
    G = zeros(1,1,Kby2);
    for i = 1:Kby2
        mat = [G_tilda(:,:,i) -G_tilda(:,:,i+Kby2-1)];
        X_tilda = null(mat);
        C1 = X_tilda(1:sz_G(2),:);
        C2 = X_tilda(sz_G(2)+1:end,:);
        null_dim = size(X_tilda, 2);
        G(1:sz_G(1),1:null_dim,i) = G_tilda(:,:,i)*C1;
        if K ~= sz_U(3)
            ii = i:K-1:sz_U(3);
            jj = i+Kby2-1:K-1:sz_U(3);
            for ind = ii
                U(1:n_r,1:null_dim,ind) = U_tilda(:,:,ind)*C1;
            end
            for ind = jj
                U(1:n_r,1:null_dim,ind) = U_tilda(:,:,ind)*C2;
            end
        else
            U(1:n_r,1:null_dim,i) = U_tilda(:,:,i)*C1;
            U(1:n_r,1:null_dim,i+Kby2-1) = U_tilda(:,:,i+Kby2-1)*C2;
        end
    end
end
% Recursive call
[G,U] = Compute_common_subspace(G, U);

% DOWNLINK MODEL
function R_S = Sum_capacity_three_cell(H, SNR, comb, K, dof)

global n_t
global n_r

%  considering only the case when M = K*N and using reciprocal channel
%  model to convert uplink equations to downlink

G = zeros(n_t, n_t, 3, 3);
W_bar = zeros(n_t, dof, 3);
V = zeros(n_t, K*dof, 3);
W = zeros(n_r, dof, 3, K);
P = zeros(K*dof, dof, 3, K);
for cell_ind = 1:3
    for other_cell_ind = 1:3
        if other_cell_ind == cell_ind, continue;end
        sz = 1;
        for user_ind = 1:K
            G(:,sz:sz+n_r-1,cell_ind,other_cell_ind) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,user_ind))';
            sz = sz + n_r;
        end
    end
end
E = (G(:,:,3,1) \ G(:,:,3,2)) * (G(:,:,1,2) \ G(:,:,1,3)) * (G(:,:,2,3) \ G(:,:,2,1));
F = G(:,:,3,2) \ G(:,:,3,1);
C = G(:,:,2,3) \ G(:,:,2,1);
[Eig_vect,~] = eig(E);
W_bar(:,:,1) = Eig_vect(:,1:dof);
W_bar(:,:,2) = F * W_bar(:,:,1);
W_bar(:,:,3) = C * W_bar(:,:,1);

for cell_ind = 1:3
    sz = 1;
    for user_ind = 1:K
        W(:,:,cell_ind,user_ind) = W_bar(sz:sz+n_r-1,:,cell_ind);
        sz = sz + n_r;
    end
end

for cell_ind = 1:3
    mat_req = zeros(n_t, 2*K*dof);
    sz = 1;
    for other_cell_ind = 1:3
        if other_cell_ind == cell_ind, continue;end
        for user_ind = 1:K
            mat_req(:,sz:sz+dof-1) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,user_ind))'*...
                W(:,:,other_cell_ind,user_ind);
            sz = sz + dof;
        end
    end
    [U,~,~] = svd(mat_req);
    sz = (2*K-1)*dof;
    V(:,:,cell_ind) = U(:,sz+1:sz+K*dof);
end

for cell_ind = 1:3
    for user_ind = 1:K
        mat_req = zeros(K*dof, dof*(K-1));
        sz = 1;
        for other_user_ind = 1:K
            if other_user_ind == user_ind, continue;end
            mat_req(:,sz:sz + dof-1) = V(:,:,cell_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,other_user_ind))'...
                *W(:,:,cell_ind,other_user_ind);
        end
        [U,~,~] = svd(mat_req);
        P(:,:,cell_ind,user_ind) = U(:,(K-1)*dof+1:end);
    end
end

R_S = 0;
Pow = 10^(SNR / 10) / K / dof;
for cell_ind = 1:3
    for user_ind = 1:K
        mat = W(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))*...
            V(:,:,cell_ind)*P(:,:,cell_ind,user_ind);
        R_S = R_S + log2(real(det(W(:,:,cell_ind,user_ind)'*W(:,:,cell_ind,user_ind) + Pow * (mat*mat')) / ...
            det(W(:,:,cell_ind,user_ind)'*W(:,:,cell_ind,user_ind))));
    end
end

% V = zeros(n_t,dof,3,K);
% for cell_ind = 1:3
%     for user_ind = 1:K
%         mat_req = zeros(n_t, (3*K-1)*dof);
%         sz = 1;
%         for other_cell_ind = 1:3
%             for other_user_ind = 1:K
%                 if other_cell_ind == cell_ind && other_user_ind == user_ind, continue;end
%                 mat_req(:,sz:sz+dof-1) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,other_user_ind))'*...
%                     W(:,:,other_cell_ind,other_user_ind);
%                 sz = sz + dof;
%             end
%         end
%         [U,~,~] = svd(mat_req);
%         sz = (2*K-1 + K - 1)*dof;
%         V(:,:,cell_ind,user_ind) = U(:,sz+1:sz+dof);
%     end
% end
% 
% R_S = 0;
% Pow = 10^(SNR / 10) / K / dof;
% for cell_ind = 1:3
%     for user_ind = 1:K
%         mat = W(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))*...
%             V(:,:,cell_ind,user_ind);
%         R_S = R_S + log2(real(det(W(:,:,cell_ind,user_ind)'*W(:,:,cell_ind,user_ind) + Pow * (mat*mat')) / ...
%             det(W(:,:,cell_ind,user_ind)'*W(:,:,cell_ind,user_ind))));
%     end
% end


% UPLINK ORIGINAL MODEL
% function R_S = Sum_capacity_three_cell(H, SNR, comb, K, dof)
% 
% global n_t
% global n_r
% 
% %  considering only the case when M = K*N
% 
% % if strcmp(dir, 'down')
% %     H = permute(H, [2,1,3,4,5]);
% %     H = conj(H);
% % end
% 
% G = zeros(n_t, n_t, 3, 3);
% W_bar = zeros(n_t, dof, 3);
% V = zeros(n_t, K*dof, 3);
% W = zeros(n_r, dof, 3, K);
% P = zeros(K*dof, dof, 3, K);
% for cell_ind = 1:3
%     for other_cell_ind = 1:3
%         if other_cell_ind == cell_ind, continue;end
%         sz = 1;
%         for user_ind = 1:K
%             G(:,sz:sz+n_r-1,cell_ind,other_cell_ind) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,user_ind));
%             sz = sz + n_r;
%         end
%     end
% end
% E = (G(:,:,3,1) \ G(:,:,3,2)) * (G(:,:,1,2) \ G(:,:,1,3)) * (G(:,:,2,3) \ G(:,:,2,1));
% F = G(:,:,3,2) \ G(:,:,3,1);
% C = G(:,:,2,3) \ G(:,:,2,1);
% [Eig_vect,~] = eig(E);
% W_bar(:,:,1) = Eig_vect(:,1:dof);
% W_bar(:,:,2) = F * W_bar(:,:,1);
% W_bar(:,:,3) = C * W_bar(:,:,1);
% 
% for cell_ind = 1:3
%     sz = 1;
%     for user_ind = 1:K
%         W(:,:,cell_ind,user_ind) = W_bar(sz:sz+n_r-1,:,cell_ind);
%         sz = sz + n_r;
%     end
% end
% 
% for cell_ind = 1:3
%     mat_req = zeros(n_t, 2*K*dof);
%     sz = 1;
%     for other_cell_ind = 1:3
%         if other_cell_ind == cell_ind, continue;end
%         for user_ind = 1:K
%             mat_req(:,sz:sz+dof-1) = H(:,:,cell_ind,other_cell_ind,comb(other_cell_ind,user_ind))...
%                 *W(:,:,other_cell_ind,user_ind);
%             sz = sz + dof;
%         end
%     end
%     [U,~,~] = svd(mat_req);
%     sz = (2*K-1)*dof;
%     V(:,:,cell_ind) = U(:,sz+1:sz+K*dof);
% end
% 
% for cell_ind = 1:3
%     for user_ind = 1:K
%         mat_req = zeros(K*dof, dof*(K-1));
%         sz = 1;
%         for other_user_ind = 1:K
%             if other_user_ind == user_ind, continue;end
%             mat_req(:,sz:sz + dof-1) = V(:,:,cell_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,other_user_ind))...
%                 *W(:,:,cell_ind,other_user_ind);
%         end
%         [U,~,~] = svd(mat_req);
%         P(:,:,cell_ind,user_ind) = U(:,(K-1)*dof+1:end);
%     end
% end
% 
% R_S = 0;
% Pow = 10^(SNR / 10);
% for cell_ind = 1:3
%     for user_ind = 1:K
%         mat = P(:,:,cell_ind,user_ind)'*V(:,:,cell_ind)'*...
%             H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))*W(:,:,cell_ind,user_ind);
%         R_S = R_S + log2(real(det(eye(dof) + Pow / K/ dof * (mat*mat'))));
%     end
% end


function R_S = Sum_capacity_cyclic_model2_b(H, SNR, comb, L, K, dof)

global n_t
global n_r

phi = zeros(n_t, K*dof, L);
U = zeros(n_r, dof, L, K);
V = zeros(K*dof, dof, L, K);

mat_temp = zeros(n_t, K*dof);
mat_temp(1:K*dof, :) = eye(K*dof);
phi(:,:,1) = mat_temp;
phi(:,:,2) = mat_temp;
if ~rem(L,2)
    i_ind = 2:2:L;
    j_ind = 1:2:L-1;
    next_i_ind = rem(i_ind + 1, L);
    next_j_ind = j_ind + 1;
else
    i_ind = 2:2:L-1;
    j_ind = 1:2:L;
    next_i_ind = i_ind + 1;
    next_j_ind = j_ind + 1;
    next_j_ind(end) = [];
end

mat_req1 = zeros(n_r*K, n_t);
mat_req2 = zeros(n_r*K, n_t);

for half_ind = 1:ceil(L/2)
    for user_ind = 1:K
        mat_req = H(:,:,i_ind(half_ind),next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind)) ...
            * phi(:,:,i_ind(half_ind));
        [U_req,~,~] = svd(mat_req);
        U(:,:,next_i_ind(half_ind),user_ind) = U_req(:,K*dof+1:(K+1)*dof);
        mat_req = H(:,:,j_ind(half_ind),next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind)) ...
            * phi(:,:,j_ind(half_ind));
        [U_req,~,~] = svd(mat_req);
        U(:,:,next_j_ind(half_ind),user_ind) = U_req(:,K*dof+1:(K+1)*dof);
    end
    
    if half_ind == ceil(L/2), break;end

    sz = 1;
    for user_ind = 1:K
        mat_req1(sz:sz+n_r-1,:) = H(:,:,i_ind(half_ind+1),next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind));
        mat_req2(sz:sz+n_r-1,:) = H(:,:,i_ind(half_ind),next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind));
        sz = sz + n_r;
    end
    phi(:,:,i_ind(half_ind+1)) = (mat_req1' / (mat_req1*mat_req1')) * mat_req2 * phi(:,:,i_ind(half_ind));
    sz = 1;
    for user_ind = 1:K
        mat_req1(sz:sz+n_r-1,:) = H(:,:,j_ind(half_ind+1),next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind));
        mat_req2(sz:sz+n_r-1,:) = H(:,:,j_ind(half_ind),next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind));
        sz = sz + n_r;
    end
    phi(:,:,j_ind(half_ind+1)) = (mat_req1' / (mat_req1*mat_req1')) * mat_req2 * phi(:,:,j_ind(half_ind));
end

% computing U_1 again due to boundary conditions
for user_ind = 1:K
    mat_req = H(:,:,2,1,comb(1,user_ind)) * phi(:,:,2);
    [U_req,~,~] = svd(mat_req);
    U(:,:,1,user_ind) = U_req(:,K*dof+1:(K+1)*dof);
end

mat =zeros(dof*(K-1), dof*K);
for cell_ind = 1:L
    for user_ind = 1:K
        other_user_set = 1:K;
        other_user_set(other_user_set == user_ind) = [];
        sz = 1;
        for other_user_ind = other_user_set
            mat(sz:sz+dof-1,:) = U(:,:,cell_ind,other_user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,other_user_ind))...
                *phi(:,:,cell_ind);
            sz = sz + dof;
        end
        X = null(mat);
        V(:,:,cell_ind,user_ind) = X(:,1:dof);
    end
end

% mat =zeros(dof*K, dof*K);
% for cell_ind = 1:L
%     sz = 1;
%     for user_ind = 1:K
%         mat(sz:sz+dof-1,:) = U(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))...
%             *phi(:,:,cell_ind);
%         sz = sz + dof;
%     end
%     mat_inv = inv(mat);
%     sz = 1;
%     for user_ind = 1:K
%         V(:,:,cell_ind,user_ind) = mat_inv(:,sz:sz+dof-1);
%         sz = sz + dof;
%     end
% end

R_S = 0;
Pow = 10^(SNR / 10) / K / dof;
for cell_ind = 1:L
    for user_ind = 1:K
        mat = U(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))...
            *phi(:,:,cell_ind)*V(:,:,cell_ind,user_ind);
        R_S = R_S + log2(real(det(eye(dof) + Pow * (mat*mat'))));
    end
end


function R_S = Sum_capacity_cyclic_model2_c(H, SNR, comb, L, K, dof)

global n_t
global n_r

phi = zeros(n_t, K*dof, L);
U = zeros(n_r, dof, L, K);
V = zeros(K*dof, dof, L, K);

mat_temp = zeros(n_t, K*dof);
mat_temp(1:K*dof, :) = eye(K*dof);
phi(:,:,1) = mat_temp;
phi(:,:,2) = mat_temp;
if ~rem(L,2)
    i_ind = 2:2:L;
    j_ind = 1:2:L-1;
    next_i_ind = rem(i_ind + 1, L);
    next_j_ind = j_ind + 1;
else
    i_ind = 2:2:L-1;
    j_ind = 1:2:L;
    next_i_ind = i_ind + 1;
    next_j_ind = j_ind + 1;
    next_j_ind(end) = [];
end

mat_req1 = zeros(dof*K, n_t);
mat_req2 = zeros(dof*K, n_t);

for half_ind = 1:ceil(L/2)
    for user_ind = 1:K
        mat_req = H(:,:,i_ind(half_ind),next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind)) ...
            * phi(:,:,i_ind(half_ind));
        [U_req,~,~] = svd(mat_req);
        U(:,:,next_i_ind(half_ind),user_ind) = U_req(:,K*dof+1:(K+1)*dof);
        mat_req = H(:,:,j_ind(half_ind),next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind)) ...
            * phi(:,:,j_ind(half_ind));
        [U_req,~,~] = svd(mat_req);
        U(:,:,next_j_ind(half_ind),user_ind) = U_req(:,K*dof+1:(K+1)*dof);
    end
    
    if half_ind == ceil(L/2), break;end
    sz = 1;
    for user_ind = 1:K
        mat_req1(sz:sz+dof-1,:) = U(:,:,next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind))'*...
            H(:,:,i_ind(half_ind+1),next_i_ind(half_ind),comb(next_i_ind(half_ind),user_ind));
        mat_req2(sz:sz+dof-1,:) = U(:,:,next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind))'*...
            H(:,:,j_ind(half_ind+1),next_j_ind(half_ind),comb(next_j_ind(half_ind),user_ind));
        sz = sz + dof;
    end
    X = null(mat_req1);
    phi(:,:,i_ind(half_ind+1)) = X(:,1:K*dof);
    X = null(mat_req2);
    phi(:,:,j_ind(half_ind+1)) = X(:,1:K*dof);
end

% computing U_1 again due to boundary conditions
for user_ind = 1:K
    mat_req = H(:,:,2,1,comb(1,user_ind)) * phi(:,:,2);
    [U_req,~,~] = svd(mat_req);
    U(:,:,1,user_ind) = U_req(:,K*dof+1:(K+1)*dof);
end

mat =zeros(dof*(K-1), dof*K);
for cell_ind = 1:L
    for user_ind = 1:K
        other_user_set = 1:K;
        other_user_set(other_user_set == user_ind) = [];
        sz = 1;
        for other_user_ind = other_user_set
            mat(sz:sz+dof-1,:) = U(:,:,cell_ind,other_user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,other_user_ind))...
                *phi(:,:,cell_ind);
            sz = sz + dof;
        end
        X = null(mat);
        V(:,:,cell_ind,user_ind) = X(:,1:dof);
    end
end

% mat =zeros(dof*K, dof*K);
% for cell_ind = 1:L
%     sz = 1;
%     for user_ind = 1:K
%         mat(sz:sz+dof-1,:) = U(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))...
%             *phi(:,:,cell_ind);
%         sz = sz + dof;
%     end
%     mat_inv = inv(mat);
%     sz = 1;
%     for user_ind = 1:K
%         V(:,:,cell_ind,user_ind) = mat_inv(:,sz:sz+dof-1);
%         sz = sz + dof;
%     end
% end

R_S = 0;
Pow = 10^(SNR / 10) / K / dof;
for cell_ind = 1:L
    for user_ind = 1:K
        mat = U(:,:,cell_ind,user_ind)'*H(:,:,cell_ind,cell_ind,comb(cell_ind,user_ind))...
            *phi(:,:,cell_ind)*V(:,:,cell_ind,user_ind);
        R_S = R_S + log2(real(det(eye(dof) + Pow * (mat*mat'))));
    end
end


function R = Compute_SINR_Capacity(H, U, V, SNR, k, l, K, L, dof)

R = 0;
% Gamma = 1;
% for j = 1:K
%     Gamma = Gamma + norm(V(:,:,j,1))^2;
% end
% Gamma = Gamma / SNR;
Pow = SNR/K/dof;
for stream_ind = 1:dof
    U_use = U(:,:,k,l);
    V_use = V(:,:,k,l);
    SINR = norm(U_use(:,stream_ind))^2;
    inter_stream_ind = 1:dof;
    inter_stream_ind(inter_stream_ind == stream_ind) = [];
    for inter_stream_ind = inter_stream_ind
        SINR = SINR + Pow*abs(U_use(:,stream_ind)'*H(:,:,k,l,l)*...
            V_use(:,inter_stream_ind))^2;
    end
    
    inter_user_ind = 1:K;
    inter_user_ind(inter_user_ind == k) = [];
    for inter_user_ind = inter_user_ind
        SINR = SINR + Pow*norm(U_use(:,stream_ind)'*H(:,:,k,l,l)*...
            V(:,:,inter_user_ind,l))^2;
    end
    
    inter_cell_ind = 1:L;
    inter_cell_ind(inter_cell_ind == l) = [];
    for inter_cell_ind = inter_cell_ind
        for user_ind = 1:K
            SINR = SINR + Pow*norm(U_use(:,stream_ind)'*H(:,:,k,l,inter_cell_ind)*...
                V(:,:,user_ind,inter_cell_ind))^2;
        end
    end
    R = R + log2(1 + Pow*abs(U_use(:,stream_ind)'*H(:,:,k,l,l)*V_use(:,stream_ind))^2 / SINR);
end


function R_F = SumCapacity_distributed_iterative(SNR_vect, BS)

global n_t
global n_r
R_F = zeros(length(SNR_vect), 1);
time = zeros(length(SNR_vect), 1);
DoF = zeros(1,BS);
fprintf('Enter the %d values of DoF Sequentially:\n', BS);
for dof_ind = 1:BS
    DoF_tmp = input('','s');
    [DoF_tmp,status] = str2num(DoF_tmp); %#ok<*ST2NM>
    if ~status
        fprintf('Wrong value input, program will terminate!\n');
        return
    end
    DoF(dof_ind) = DoF_tmp;
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
        nusers = 1;
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            R_final(avg_sim_ind) = Sum_capacity_distr_iter(H, SNR_vect(SNR_ind), DoF);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation at SNR = %2.1f dB done, time taken = %f\n', SNR_vect(SNR_ind), time(SNR_ind));
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_Test(SNR_vect, L, users_num)

global n_t
global n_r
R_F = zeros(length(SNR_vect), 1);
time = zeros(length(SNR_vect), 1);
DoF = zeros(1,L);
% K = 2;
fprintf('Enter the %d values of DoF Sequentially:\n', L);
for dof_ind = 1:L
    DoF_tmp = input('','s');
    [DoF_tmp,status] = str2num(DoF_tmp); %#ok<*ST2NM>
    if ~status
        fprintf('Wrong value input, program will terminate!\n');
        return
    end
    DoF(dof_ind) = DoF_tmp;
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
        nusers = users_num;
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            H = (randn(n_r,n_t,nusers,L,L) + 1i*randn(n_r,n_t,nusers,L,L)) / sqrt(2);
            R_final(avg_sim_ind) = Sum_capacity_distr_iter_GG(H, SNR_vect(SNR_ind), DoF, nusers);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation at SNR = %2.1f dB done, time taken = %f\n', SNR_vect(SNR_ind), time(SNR_ind));
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_GG(SNR_vect, L, users_num)

global n_t
global n_r

R_F = zeros(length(SNR_vect), length(users_num));
time = zeros(length(SNR_vect), length(users_num));
fprintf('Enter the value of DoF for each user in each cell :\n');
DoF_str = input('','s');
[DoF,status] = str2num(DoF_str);
if ~status
    fprintf('Wrong value input, program will terminate!\n');
    return
end
K = min(floor((n_t/DoF-1)/(L-1)), 1+floor((n_r/DoF-1)/(L-1)));
if K<=1
    fprintf('Unsupported value input (K<=1), program will terminate!\n');
    return
end 
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
    for mi = 1:length(users_num)
        nusers = users_num(mi);
        R_final = zeros(N_sim,1);    
%         str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,L);
%         eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

%             H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            H = (randn(n_r,n_t,nusers,L,L) + sqrt(-1)*randn(n_r,n_t,nusers,L,L))/sqrt(2);
            
            user_comb = zeros(L,K);
            for cell_ind = 1:L
                user_comb(cell_ind,:) = getInitialUserSet('old',H(:,:,:,cell_ind,cell_ind), K, 1:nusers);
            end

            R_temp = 0;
            for cell_ind = 1:L
                Ent_temp = 0;
                user_set = 1:nusers;
                for k = 1:K
                    T_set = user_set;
                    other_user_ind = 1:K;
                    other_user_ind(other_user_ind == k) = [];
                    other_users = user_comb(cell_ind,other_user_ind);
                    T_set(T_set == other_users) = [];
                    comb_temp = user_comb;
                    for user_ind = T_set
                        comb_temp(cell_ind,k) = user_ind;
                        Ent = compute_selection_metric(H, comb_temp, K, L, DoF, [k,cell_ind]);
                        if Ent > Ent_temp
                            Ent_temp = Ent;
                            maxValInd = user_ind;
                        end
                    end
                    comb_temp = user_comb;
                    comb_temp(cell_ind,k) = maxValInd;
                    R_S = Sum_capacity_grouping(H, SNR_vect(SNR_ind), comb_temp, K,L, DoF);
                    if R_S > R_temp
                        user_comb = comb_temp;
                        R_temp = R_S;
                    end
                end
            end
            R_final(avg_sim_ind) = R_temp;
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind,mi) = mean(R_final);
        time(SNR_ind,mi) = toc;
        fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', nusers, SNR_vect(SNR_ind), time(SNR_ind,mi));
    end
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function met_val = compute_selection_metric(H, comb, K, L, dof, user_id)

global n_t
global n_r

U = zeros(n_r,dof,K,L);
G = zeros(n_t,dof,L);

for cell_ind = 1:L
    next_cell_ind = cell_ind + 1;
    if cell_ind + 1 > L, next_cell_ind = 1;end
    if K <= 3
        mat_req = zeros(K*n_t,n_t + K*n_r);
        sz_eye = 1;
        sz_H_colm = n_t;
        sz_H_row = 1;
        for i = 1:K
            mat_req(sz_eye:sz_eye+n_t-1,1:n_t) = eye(n_t);
            mat_req(sz_H_row:sz_H_row+n_t-1,sz_H_colm+1:sz_H_colm+n_r) = -H(:,:,comb(next_cell_ind,i),next_cell_ind,cell_ind)';
            sz_eye = sz_eye + n_t;
            sz_H_row = sz_H_row + n_t;
            sz_H_colm = sz_H_colm + n_r;
        end
        X = null(mat_req);
        G(:,:,cell_ind) = X(1:n_t,1:dof);
        sz = n_t;
        for i = 1:K
            U(:,:,i,next_cell_ind) = X(sz+1:sz+n_r,1:dof);
            sz = sz + n_r;
        end
    else
        G_tilda = zeros(n_t,n_r,K);
        U_tilda = zeros(n_r,n_r,K);
        for k = 1:K
            mat = [eye(n_t) -H(:,:,comb(next_cell_ind,k),next_cell_ind,cell_ind)'];
            X_tilda = null(mat);
            G_tilda(:,:,k) = X_tilda(1:n_t,:);
            U_tilda(:,:,k) = X_tilda(n_t+1:end,:);
        end
        [G_get,U_get] = Compute_common_subspace(G_tilda, U_tilda);
        G(:,:,cell_ind) = G_get;
        U(:,:,:,next_cell_ind) = U_get;
    end
end

next_cell_ind = user_id(2) + 1;
if user_id(2) + 1 > L, next_cell_ind = 1;end
H_A = H(:,:,comb(user_id(2),user_id(1)),user_id(2),user_id(2))' * U(:,:,user_id(1),user_id(2));
cell_set = 1:L;
cell_set(cell_set == next_cell_ind) = [];
H_B_hor = zeros(n_t,(K*(L-1)-1)*dof);
H_B_ver = zeros((K*(L-1)-1)*n_t,dof);
sz_hor = 1;
sz_ver = 1;
for cell_ind = cell_set
    for user_ind = 1:K
        if user_ind == user_id(1) && cell_ind == user_id(2),continue;end
        H_B_hor(:,sz_hor:sz_hor+dof-1) = H(:,:,comb(cell_ind,user_ind),cell_ind,user_id(2))'*...
            U(:,:,user_ind,cell_ind);
        H_B_ver(sz_ver:sz_ver+n_t-1,:) = H(:,:,comb(cell_ind,user_ind),cell_ind,user_id(2))'*...
            U(:,:,user_ind,cell_ind);
        sz_hor = sz_hor + dof;
        sz_ver = sz_ver + n_t;
    end
end

H_B_hor = [H_B_hor,G(:,:,user_id(2))];
H_B_ver = [H_B_ver;G(:,:,user_id(2))];
met_val = cho_dist(H_A,H_B_hor);
% P = 10;
% Ent1 = log2(real(det(eye(dof) + P/n_r*(H_A'*H_A))));
% Ent2 = log2(real(det(eye(dof) + P/n_r*(H_B_ver'*H_B_ver))));
% Ent12 = log2(real(det(eye(dof) + P/n_r*(H_A'*H_A) + P/n_r*(H_B_ver'*H_B_ver))));
% met_val = 2*Ent12 - Ent1 - Ent2;


function R_F = SumCapacity_iterative_max_SINR(SNR_vect, BS)

global n_t
global n_r
R_F = zeros(length(SNR_vect), 1);
time = zeros(length(SNR_vect), 1);
DoF = zeros(1,BS);
fprintf('Enter the %d values of DoF Sequentially:\n', BS);
for dof_ind = 1:BS
    DoF_tmp = input('','s');
    [DoF_tmp,status] = str2num(DoF_tmp); %#ok<*ST2NM>
    if ~status
        fprintf('Wrong value input, program will terminate!\n');
        return
    end
    DoF(dof_ind) = DoF_tmp;
end
N_sim = 1000;
dlg = ProgressDialog();
for SNR_ind = 1:length(SNR_vect)
        nusers = 1;
        R_final = zeros(N_sim,1);    
        str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
        eval(['load ' str]);
        dlg.FractionComplete = 0;
        dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
        tic
        for avg_sim_ind = 1:N_sim

            H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
            R_final(avg_sim_ind) = Sum_capacity_iter_mSINR(H, SNR_vect(SNR_ind), DoF);
            if ~rem(avg_sim_ind,10)
                dlg.FractionComplete = avg_sim_ind / N_sim;
                dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
            end
        end 
        R_F(SNR_ind) = mean(R_final);
        time(SNR_ind) = toc;
        fprintf('Simulation at SNR = %2.1f dB done, time taken = %f\n', SNR_vect(SNR_ind), time(SNR_ind));
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_S = Sum_capacity_distr_iter(H, SNR, dof)

global n_t
global n_r

R_S = 0;
SNR = 10 ^ (SNR / 10);
K = length(dof);
V = cell(1,K);
U = cell(1,K);
for i = 1:K
    V{i} = zeros(n_t,dof(i));
    V{i}(1:dof(i),:) = eye(dof(i));
    U{i} = zeros(dof(i),n_r);
    U{i}(:,1:dof(i)) = eye(dof(i));
end
% Inter_leak_prev = 0;
niter = 0;
while(1)
    for k = 1:K
        Q = zeros(n_r,n_r);
        for j = 1:K
            if j == k, continue;end
            Q = Q + SNR/dof(j) * H(:,:,1,k,j)*V{j}*V{j}'*H(:,:,1,k,j)';
        end
        [eig_vect,eig_val] = eig(Q);
        [~, ind_sort] = sort(diag(eig_val),'ascend');
        U{k} = eig_vect(:,ind_sort(1:dof(k)));
    end
    
    for k = 1:K
        Q = zeros(n_t,n_t);
        for j = 1:K
            if j == k, continue;end
            Q = Q + SNR / dof(j) * H(:,:,1,j,k)'*U{j}*U{j}'*H(:,:,1,j,k);
        end
        [eig_vect,eig_val] = eig(Q);
        [~, ind_sort] = sort(diag(eig_val),'ascend');
        V{k} = eig_vect(:,ind_sort(1:dof(k)));
    end
    Inter_leak = 0;
    for k = 1:K
        Q = zeros(n_r,n_r);
        for j = 1:K
            if j == k, continue;end
            Q = Q + SNR/dof(j) * H(:,:,1,k,j)*V{j}*V{j}'*H(:,:,1,k,j)';
        end
%         [~,eig_val] = eig(Q);
%         [eig_val] = sort(diag(eig_val),'ascend');
%         sum(eig_val(1:dof(k)))
%         Q = zeros(n_r,n_r);
        Inter_leak = Inter_leak + trace(U{k}'*Q*U{k});
    end
    Inter_leak = real(Inter_leak);
    if abs(Inter_leak) < 1e-4
        break;
    end
    niter = niter + 1;
%     if niter > 50,
%         disp(Inter_leak);
%         break;
%     end
end
for k = 1:K
    H_k = U{k}'*H(:,:,1,k,k)*V{k};
    R_S = R_S + log2(real(det(eye(dof(k)) + SNR/dof(k)*(H_k*H_k'))));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Below formulation is also correct and is giving R / log(P) -> 3
% for i = 1:K
%     for j = 1:dof
%         SINR = 1;
%         for k = 1:dof
%             if k == j, continue;end
%             SINR = SINR + SNR/dof*abs(U{i}(:,j)'*H(:,:,1,i,i)*V{i}(:,k))^2;
%         end
%         R_S = R_S + log2(1 + SNR/dof(i)*abs(U{i}(:,j)'*H(:,:,1,i,i)*V{i}(:,j))^2 / SINR);
%     end
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function R_S = Sum_capacity_distr_iter_GG(H, SNR, dof, K)

global n_t
global n_r

R_S = 0;
SNR = 10 ^ (SNR / 10);
L = length(dof);
V = cell(K,L);
U = cell(K,L);
for i = 1:L
    for j = 1:K
        V{j,i} = zeros(n_t,dof(i));
        V{j,i}(1:dof(i),:) = eye(dof(i));
        U{j,i} = zeros(dof(i),n_r);
        U{j,i}(:,1:dof(i)) = eye(dof(i));
    end
end
% Inter_leak_prev = 0;
niter = 0;
inter_len = 500;
Inter_leak = zeros(inter_len,1);
while(1)
    for k = 1:L
        for u_i = 1:K
            Q = zeros(n_r,n_r);
            for j = 1:L
                for u_ij = 1:K
                    if j == k && u_ij == u_i, continue;end
                    Q = Q + SNR/dof(j) * H(:,:,u_i,k,j)*V{u_ij,j}*V{u_ij,j}'*H(:,:,u_i,k,j)';
                end
            end
            [eig_vect,eig_val] = eig(Q);
            [~, ind_sort] = sort(diag(eig_val),'ascend');
            U{u_i,k} = eig_vect(:,ind_sort(1:dof(k)));
        end
    end
    
    for k = 1:L
        for u_i = 1:K
            Q = zeros(n_t,n_t);
            for j = 1:L
                for u_ij = 1:K
                    if j == k && u_ij == u_i, continue;end
                    Q = Q + SNR / dof(j) * H(:,:,u_ij,j,k)'*U{u_ij,j}*U{u_ij,j}'*H(:,:,u_ij,j,k);
                end
            end
            [eig_vect,eig_val] = eig(Q);
            [~, ind_sort] = sort(diag(eig_val),'ascend');
            V{u_i,k} = eig_vect(:,ind_sort(1:dof(k)));
        end
    end
    Inter_leak_temp = 0;
    for k = 1:L
        for u_i = 1:K
            Q = zeros(n_r,n_r);
            for j = 1:L
                for u_ij = 1:K
                    if j == k && u_ij == u_i, continue;end
                    Q = Q + SNR/dof(j) * H(:,:,u_i,k,j)*V{u_ij,j}*V{u_ij,j}'*H(:,:,u_i,k,j)';
                end
            end
            Inter_leak_temp = Inter_leak_temp + trace(U{u_i,k}'*Q*U{u_i,k});
        end
    end
    Inter_leak(niter+1) = real(Inter_leak_temp);
%     if abs(Inter_leak(niter+1)) < 1e-4 || niter == inter_len-1
    if  niter == inter_len-1
        figure;
        plot(1:inter_len,Inter_leak);grid;
        break;
    end
    niter = niter + 1;
end
for k = 1:L
    H_k = U{k}'*H(:,:,1,k,k)*V{k};
    R_S = R_S + log2(real(det(eye(dof(k)) + SNR/dof(k)*(H_k*H_k'))));
end


function R_S = Sum_capacity_iter_mSINR(H, SNR, dof)

global n_t
global n_r

SNR = 10 ^ (SNR / 10);
K = length(dof);
V = cell(1,K);
U = cell(1,K);
for i = 1:K
    V{i} = zeros(n_t,dof(i));
    V{i}(1:dof(i),:) = eye(dof(i));
    U{i} = zeros(n_r,dof(i));
end
% U_prev = U;
% V_prev = V;
niter = 0;
SINR_prev = zeros(K, max(dof));
while(1)
    for k = 1:K
        for l = 1:dof(k)
            B_kl = eye(n_r);
            for j = 1:K
                mat_use = zeros(n_r,n_r);
                for d = 1:dof(j)
                    mat_use = mat_use + H(:,:,1,k,j)*V{j}(:,d)*V{j}(:,d)'*...
                        H(:,:,1,k,j)';
                end
                B_kl = B_kl + SNR/dof(j) * mat_use;
            end
            B_kl = B_kl - SNR/dof(k)*H(:,:,1,k,k)*V{k}(:,l)*V{k}(:,l)'*...
                H(:,:,1,k,k)';
            vect_use = (B_kl\ H(:,:,1,k,k))*V{k}(:,l);
            U{k}(:,l) = vect_use / norm(vect_use);
        end
    end
    
    for k = 1:K
        for l = 1:dof(k)
            B_kl = eye(n_r);
            for j = 1:K
                mat_use = zeros(n_t,n_t);
                for d = 1:dof(j)
                    mat_use = mat_use + H(:,:,1,j,k)'*U{j}(:,d)*U{j}(:,d)'*...
                        H(:,:,1,j,k);
                end
                B_kl = B_kl + SNR/dof(j) * mat_use;
            end
            B_kl = B_kl - SNR/dof(k)*H(:,:,1,k,k)'*U{k}(:,l)*U{k}(:,l)'*...
                H(:,:,1,k,k);
            vect_use = (B_kl\ (H(:,:,1,k,k)'))*U{k}(:,l);
            V{k}(:,l) = vect_use / norm(vect_use);
        end
    end
    
%     iter_norm = 0;
%     for i = 1:K
%         mat = U_prev{k} - U{k};
%         iter_norm = iter_norm + norm(mat,'fro')^2;
%         mat = V_prev{k} - V{k};
%         iter_norm = iter_norm + norm(mat,'fro')^2;
%     end
%     if iter_norm < 1e-7,
%         break;
%     end
    SINR = zeros(K,max(dof));
    for k = 1:K
        for l = 1:dof(k)
            B_kl = eye(n_r);
            for j = 1:K
                mat_use = zeros(n_r,n_r);
                for d = 1:dof(j)
                    mat_use = mat_use + H(:,:,1,k,j)*V{j}(:,d)*V{j}(:,d)'*...
                        H(:,:,1,k,j)';
                end
                B_kl = B_kl + SNR/dof(j) * mat_use;
            end
            B_kl = B_kl - SNR/dof(k)*H(:,:,1,k,k)*V{k}(:,l)*V{k}(:,l)'*...
                H(:,:,1,k,k)';
            SINR(k,l) = SNR / dof(k) * abs(U{k}(:,l)'*H(:,:,1,k,k)*V{k}(:,l))^2;
            SINR(k,l) = SINR(k,l) / (U{k}(:,l)'*B_kl*U{k}(:,l));
        end
    end
    SINR = real(SINR);
    if norm(SINR - SINR_prev, 'fro')^2 < 1e-7
        break;
    end
    
    niter = niter + 1;
    if niter > 100
        break;
    end
%     U_prev = U;
%     V_prev = V;
    SINR_prev = SINR;
end

R_S = sum(sum(log2(1 + SINR)));

% for k = 1:K
%     for l = 1:dof(k)
%         B_kl = eye(n_r);
%         for j = 1:K
%             mat_use = zeros(n_r,n_r);
%             for d = 1:dof(j)
%                 mat_use = mat_use + H(:,:,1,k,j)*V{j}(:,d)*V{j}(:,d)'*...
%                     H(:,:,1,k,j)';
%             end
%             B_kl = B_kl + SNR/dof(j) * mat_use;
%         end
%         B_kl = B_kl - SNR/dof(k)*H(:,:,1,k,k)*V{k}(:,l)*V{k}(:,l)'*...
%             H(:,:,1,k,k)';
%         SINR = SNR / dof(k) * abs(U{k}(:,l)'*H(:,:,1,k,k)*V{k}(:,l))^2;
%         SINR = SINR / (U{k}(:,l)'*B_kl*U{k}(:,l));
%         R_S = R_S + log2(1 + real(SINR));
%     end
% end


function R_S = sum_capacity_R_sigma(H, SNR, K_set)

% global n_t
R_S = 0;
No_BS = length(K_set);
if No_BS ~= 3
    disp('The post pre processing design restricted to only 3 BS');
    return;
end

% [U, V] = compute_post_pre_processing_matrices(H, K_set, SNR);
% [U, V] = compute_post_pre_processing_matrices_MMSE_opt(H, K_set, SNR);
R_S = get_sum_capacity_ZF_optimal(H, K_set, SNR);

% for cell_ind = 1:No_BS
%     for stream_ind = 1:n_t/2
%         interf = 0;
%         for inter_cell_ind = 1:No_BS
%             if inter_cell_ind == cell_ind,continue;end
%             interf = interf + norm(U{cell_ind}(stream_ind,:)*...
%                 H(:,:,K_set(cell_ind),cell_ind,inter_cell_ind)*V{inter_cell_ind})^2;
%         end
%         for inter_stream_ind = 1:n_t/2
%             if inter_stream_ind == stream_ind,continue;end
%             interf = interf + abs(U{cell_ind}(stream_ind,:)*...
%                 H(:,:,K_set(cell_ind),cell_ind,cell_ind)*V{cell_ind}(:,inter_stream_ind))^2;
%         end
%         SINR = abs(U{cell_ind}(stream_ind,:)*H(:,:,K_set(cell_ind),cell_ind,cell_ind)*V{cell_ind}(:,stream_ind))^2 / ...
%             (norm(U{cell_ind}(stream_ind,:))^2 + interf);
%         R_S = R_S + log2(1 + SINR);
%     end
% end
% for i = 1:3
%     mat = eye(n_t/2);
%     for j = 1:3
%         if j == i, continue;end
%         mat = mat + U{i}*H(:,:,1,i,j)*V{j}*V{j}'*H(:,:,1,i,j)'*(U{i}');
%     end
%     sum_mat = eye(n_t/2) + (U{i}*H(:,:,1,i,i)*V{i}*V{i}'*H(:,:,1,i,i)'*(U{i}')) / mat;
%     R_S = R_S + log2(real(det(sum_mat)));
% end

% for i = 1:3
%     mat = eye(n_t);
%     for j = 1:3
%         if j == i, continue;end
%         mat = mat + H(:,:,1,i,j)*V{j}*V{j}'*H(:,:,1,i,j)';
%     end
%     sum_mat = eye(n_t) + (H(:,:,1,i,i)*V{i}*V{i}'*H(:,:,1,i,i)') / mat;
%     R_S = R_S + log2(real(det(sum_mat)));
% end


function [U, V] = compute_post_pre_processing_matrices(H, K_set, SNR)

global n_t

SNR = 10 ^ (SNR / 10);
No_BS = length(K_set);
V = cell(No_BS, 1);
U = cell(No_BS,1);
E = (H(:,:,K_set(3),3,1)\H(:,:,K_set(3),3,2))*(H(:,:,K_set(1),1,2)\...
    H(:,:,K_set(1),1,3))*(H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1));
[eig_vects,~] = eig(E);
V{1} = eig_vects(:,1:n_t/2);
V{2} = (H(:,:,K_set(3),3,2)\H(:,:,K_set(3),3,1))*V{1};
V{3} = (H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1))*V{1};

V{1} = V{1}* sqrt(SNR  / trace(V{1}*V{1}'));
V{2} = V{2}* sqrt(SNR  / trace(V{2}*V{2}'));
V{3} = V{3}* sqrt(SNR  / trace(V{3}*V{3}'));

% V{1} = V{1}* sqrt(SNR / 3 / trace(V{1}*V{1}'));
% V{2} = V{2}* sqrt(SNR / 3 / trace(V{2}*V{2}'));
% V{3} = V{3}* sqrt(SNR / 3 / trace(V{3}*V{3}'));

for cell_ind = 1:No_BS
    H_e = [H(:,:,K_set(cell_ind),cell_ind,1)*V{1}, ...
        H(:,:,K_set(cell_ind),cell_ind,2)*V{2},...
        H(:,:,K_set(cell_ind),cell_ind,3)*V{3}];
%     H_e = [H(:,:,K_set(cell_ind),cell_ind,1)*V{1}; ...
%         H(:,:,K_set(cell_ind),cell_ind,2)*V{2};...
%         H(:,:,K_set(cell_ind),cell_ind,3)*V{3}];
    U_bar = (H_e'*H_e + eye(3*n_t/2))\(H_e');
    U{cell_ind} = U_bar((cell_ind-1)*n_t/2+1:cell_ind*n_t/2,:);
%     U{cell_ind} = U_bar(:,(cell_ind-1)*n_t+1:cell_ind*n_t);
end


function R_F = SumCapacity_conventional(SNR_vect, BS, users_num)

global n_t
global n_r
if users_num ~= 1, return;end
% SNR_vect = 10.^(SNR_vect/10);
R_F = zeros(length(SNR_vect), 1);
time = zeros(length(SNR_vect), 1);
N_sim = 1000;
dlg = ProgressDialog();
nusers = 1;
str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
eval(['load ' str]);
for SNR_ind = 1:length(SNR_vect)
    dlg.FractionComplete = 0;
    dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f, %d%% done', nusers,SNR_vect(SNR_ind),0);
    R_final = zeros(N_sim,1);    
    tic
    for avg_sim_ind = 1:N_sim

        H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
%         R_final(avg_sim_ind) = sum_capacity_conventional(H, [1,1,1], SNR_vect(SNR_ind));
        R_final(avg_sim_ind) = sum_capacity_R_sigma(H, SNR_vect(SNR_ind), [1,1,1]);
        if ~rem(avg_sim_ind,10)
            dlg.FractionComplete = avg_sim_ind / N_sim;
            dlg.StatusMessage = sprintf('Number of users = %d, SNR = %2.1f dB, %d%% done',nusers,SNR_vect(SNR_ind),avg_sim_ind*100/N_sim);
        end
    end
    R_F(SNR_ind,1) = mean(R_final);
    time(SNR_ind,1) = toc;
    fprintf('Simulation for no. of users = %d at SNR = %2.1f dB done, time taken = %f\n', users_num(1), SNR_vect(SNR_ind), time(SNR_ind,1));
end
if numel(time) > 1
    fprintf('Total time taken = %f\n', sum(sum(time)));
end


function R_F = SumCapacity_TDMA(SNR, BS, users_num)

global n_t
global n_r
R_F = zeros(1, length(users_num));
time = zeros(1, length(users_num));
N_sim = 1000;
dlg = ProgressDialog();
for mi = 1:length(users_num)
    nusers = users_num(mi);
    R_final = zeros(N_sim,1);    
    str = sprintf('Channel_Mat/H_all_%d_%d_%d_%d',n_t,n_r,nusers,BS);
    eval(['load ' str]);
    dlg.FractionComplete = 0;
    dlg.StatusMessage = sprintf('Number of users = %d, %d%% done', nusers,0);
    tic
    for avg_sim_ind = 1:N_sim
        
        H = H_all(:,:,:,:,:,avg_sim_ind); %#ok<*NODEF>
        R_temp = 0;
        for i = 1:BS
            for j = 1:nusers
                R_S = sum_capacity_R_TDMA(H, SNR, [i,j]);
                if R_S > R_temp
                    R_temp = R_S;
                end
            end
        end
        R_final(avg_sim_ind) = R_temp;
        if ~rem(avg_sim_ind,10)
            dlg.FractionComplete = avg_sim_ind / N_sim;
            dlg.StatusMessage = sprintf('Number of users = %d, %d%% done',nusers,avg_sim_ind*100/N_sim);
        end
    end
    R_F(mi) = mean(R_final);
    time(mi) = toc;
    fprintf('Simulation for no. of users = %d done, time taken = %f\n', users_num(mi), time(mi));
end
if length(time) > 1
    fprintf('Total time taken = %f\n', sum(time));
end


function [U, V] = compute_post_pre_processing_matrices_MMSE_opt(H, K_set, SNR)

global n_t
SNR = 10^(SNR/10);
No_BS = length(K_set);
V = cell(No_BS, 1);
U = cell(No_BS, 1);
E = (H(:,:,K_set(3),3,1)\H(:,:,K_set(3),3,2))*(H(:,:,K_set(1),1,2)\...
    H(:,:,K_set(1),1,3))*(H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1));
[eig_vects,~] = eig(E);
V{1} = eig_vects(:,1:n_t/2);
V{2} = (H(:,:,K_set(3),3,2)\H(:,:,K_set(3),3,1))*V{1};
V{3} = (H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1))*V{1};

P = ones(1,3)*SNR;
Q = cell(1,No_BS);
[Q{1},~] = qr(V{1},0);
[Q{2},~] = qr(V{2},0);
[Q{3},~] = qr(V{3},0);
M_bar = cell(1,No_BS);
for i = 1:No_BS
    mat = eye(n_t)*n_t /(2*P(i)) ;
    for j = 1:No_BS
        if j == 1, continue;end
        mat = mat + H(:,:,K_set(i),i,j)*Q{j}*Q{j}'*H(:,:,K_set(i),i,j)' * P(j) / P(i);
    end
    M_bar{i} = (H(:,:,K_set(i),i,i)*Q{i})' * mat^-1;
end
for i = 1:No_BS
    for j = 1:No_BS
        mat = eye(n_t);
        if j == i, continue;end
        mat = mat + H(:,:,K_set(i),i,j)*Q{j}*Q{j}'*H(:,:,K_set(i),i,j)' * P(j);
    end
    A = M_bar{i}*mat*M_bar{i}';
    A(1:n_t/2+1:end) = real(A(1:n_t/2+1:end));
    L = chol(A, 'lower');
    H_ms = L \ M_bar{i} * H(:,:,K_set(i),i,i)*Q{i};
    [U_gg,~,V_gg] = svd(H_ms);
    U{i} = U_gg' * (L \ M_bar{i});
    V{i} = sqrt(2*P(i)/n_t) * Q{i} * V_gg;
end

function R_S = get_sum_capacity_ZF_optimal(H, K_set, SNR)

global n_t
R_S = 0;
SNR = 10^(SNR/10);
No_BS = length(K_set);
V = cell(No_BS, 1);
E = (H(:,:,K_set(3),3,1)\H(:,:,K_set(3),3,2))*(H(:,:,K_set(1),1,2)\...
    H(:,:,K_set(1),1,3))*(H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1));
[eig_vects,~] = eig(E);
V{1} = eig_vects(:,1:n_t/2);
V{2} = (H(:,:,K_set(3),3,2)\H(:,:,K_set(3),3,1))*V{1};
V{3} = (H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1))*V{1};

P = ones(1,3)*SNR;
Q = cell(1,No_BS);
[Q{1},~] = qr(V{1},0);
[Q{2},~] = qr(V{2},0);
[Q{3},~] = qr(V{3},0);
M_bar = cell(1,No_BS);
for i = 1:No_BS
    for j = 1:No_BS
        if j == i, continue;end
        [U_i,~,~] = svd(H(:,:,K_set(i),i,j)*Q{j});
        break;
    end
    M_bar{i} = U_i(:,n_t/2+1:end)';
end
for i = 1:No_BS
    H_zf_i = M_bar{i}*H(:,:,K_set(i),i,i)*Q{i};
    [~,L_i,V_i] = svd(H_zf_i);
    lambda = diag(L_i).^2;
    pl = Water_filling(lambda, P(i));
    sigma = diag(pl);
    C_zf = V_i * (sigma^0.5);
    R_S = R_S + log2(real(det(eye(n_t/2) + L_i*V_i'*(C_zf*C_zf')*V_i*L_i')));
end


function R_F = sum_capacity_conventional(H, K_set, SNR)

global n_t

R_F = 0;
V = cell(3,1);
E = (H(:,:,K_set(3),3,1)\H(:,:,K_set(3),3,2))*(H(:,:,K_set(1),1,2)\...
    H(:,:,K_set(1),1,3))*(H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1));
[eig_vects,~] = eig(E);
V{1} = eig_vects(:,1:n_t/2);
V{2} = (H(:,:,K_set(3),3,2)\H(:,:,K_set(3),3,1))*V{1};
V{3} = (H(:,:,K_set(2),2,3)\H(:,:,K_set(2),2,1))*V{1};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% equal power division among all BS (number = 3)

% V{1} = V{1}* sqrt(SNR / 3 / trace(V{1}*V{1}'));
% V{2} = V{2}* sqrt(SNR / 3 / trace(V{2}*V{2}'));
% V{3} = V{3}* sqrt(SNR / 3 / trace(V{3}*V{3}'));

% for i = 1:3
%     mat = eye(n_t);
%     for j = 1:3
%         if j == i, continue;end
%         mat = mat + H(:,:,1,i,j)*V{j}*V{j}'*H(:,:,1,i,j)';
%     end
%     sum_mat = eye(n_t) + (H(:,:,1,i,i)*V{i}*V{i}'*H(:,:,1,i,i)') / mat;
%     R_F = R_F + log2(real(det(sum_mat)));
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% For paper by Cadambe, theoretical Inter_Align
V{1} = GSO(V{1});
V{2} = GSO(V{2});
V{3} = GSO(V{3});
U = cell(1,3);
for i = 1:3
    for j = 1:3
        if j == i, continue;end
        [Us,~,~] = svd(H(:,:,1,i,j)*V{j});
        U{i} = Us(:,n_t/2+1:end);
        break
    end
end

for k = 1:3
    H_kk = U{k}'*H(:,:,1,k,k)*V{k};
    R_F = R_F + log2(real(det(eye(n_t/2) + SNR*2/n_t * (H_kk*H_kk'))));
end


function R_S = sum_capacity_R_TDMA(H, SNR, User_ID)

H_req = H(:,:,User_ID(2),User_ID(1),User_ID(1));
lambda = eig(H_req*H_req');
pl = Water_filling(lambda, SNR);
R_S = sum(log2(1 + lambda.*pl));


function [pl] = Water_filling(lambda, P)

[lambda idx] = sort(lambda, 'descend');
lambda = lambda(find(lambda > 1e-7));      %#ok<FNDSB> % ignoring non-positive eigenvalues
pl = -1;
try
    while (min(pl) < 0)
        mu = (P + sum(1 ./ lambda)) / length(lambda);
        pl = mu - 1 ./ lambda;
        lambda = lambda(1:end-1);
    end
catch %#ok<*CTCH>
    disp('There exists no water filling level for the input eigenvalues. Check your data and try again');
end
pl = [pl; zeros(length(idx) - length(pl), 1)]; % assigning zero power for weak eigen-modes
pl(idx) = pl; % rearranging the power levels


function d = cho_dist(A,B)

% A = GSO(A.');
% A = A.';
% B = GSO(B.');
% B = B.';
A = GSO(A);
B = GSO(B);
d = norm((A*A') - (B*B'), 'fro');
% d = norm((A'*A) - (B'*B), 'fro');


function[Q] =  GSO(A)

[~,n] = size(A);
% compute Q R using Gram-Schmidt
for j = 1:n
   v = A(:,j);
   for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
   end
   R(j,j) = norm(v);
   Q(:,j) = v/R(j,j);
end


function A = allcomb_gg(comb, ni)

ii = ni:-1:1 ;
[A{ii}] = ndgrid(comb) ;   
A = reshape(cat(ni+1,A{:}),[],ni) ;