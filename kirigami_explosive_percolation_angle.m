% Simulation for explosive rigidity percolation using the angle-based method
% 
% - Rigidify a quad or hole by fixing its angle structure each step 
%   based on some selection rule
% - Check the angle constraints to update the rigidity of other tiles/holes
% - Count the total number of rigid tiles and holes
%
% Reference:
% G. P. T. Choi, L. Liu, L. Mahadevan, "Explosive rigidity perolcation in
% kirigami", preprint, 2022.
% 
% Copyright (c) 2022,  Gary P. T. Choi, L. Liu, L. Mahadevan

mkdir('result_angle');
addpath('result_angle');

%% Setup

% kirigami pattern size
L_all = [5,10,15,20];

% number of choices in each step
k_all = [1,2,5,10,25];

% selection rule (1: most efficient, 2: least efficient)
rule_all = [1,2]; 

% number of simulations for each setup
n_sim = 200;

%% Run the simulations (may take several hours for large L and large k)

for rule = rule_all
    for L = L_all
        
        N_max = L^2+(L-1)^2;
        N_rigid_all = zeros(N_max,n_sim);
            
        for k = k_all
            parfor jjj = 1:n_sim
                %% start a simulation
                
                A = zeros(L,L); % tile rigidity
                B = zeros(L-1,L-1); % hole rigidity
                candidates = 1:N_max;
                chosen = [];
                
                for iii = 1:N_max
                    if sum(sum(A))+sum(sum(B)) ~= N_max
                        % k choices
                        if length(candidates)~=1
                            ids_batch = randsample(candidates,min(length(candidates),k));
                        else
                            ids_batch = candidates;
                        end
                        if k > 1
                            %% consider the outcome for each choice
                            N_rigid_trial = zeros(min(length(candidates),k),1);
                            for kkk = 1:min(length(candidates),k)
                                A_trial = A;
                                B_trial = B;
                                if ids_batch(kkk) <= L^2
                                    A_trial(ids_batch(kkk)) = 1;
                                else
                                    B_trial(ids_batch(kkk)-L^2) = 1;
                                end

                                A_temp = zeros(L+2,L+2);
                                B_temp = zeros(L+1,L+1);
                                A_temp(2:end-1,2:end-1) = A_trial;
                                B_temp(2:end-1,2:end-1) = B_trial;

                                % check if a hole becomes rigid
                                for i = 1:L-1
                                    for j = 1:L-1
                                        B_trial(i,j) = max([B_temp(i+1,j+1),...
                                            B_temp(i,j+1)*A_temp(i+1,j+1)*A_temp(i+1,j+1+1),...
                                            B_temp(i+1+1,j+1)*A_temp(i+1+1,j+1)*A_temp(i+1+1,j+1+1),...
                                            B_temp(i+1,j-1+1)*A_temp(i+1,j+1)*A_temp(i+1+1,j+1),...
                                            B_temp(i+1,j+1+1)*A_temp(i+1,j+1+1)*A_temp(i+1+1,j+1+1)]);
                                    end
                                end

                                % check if a quad becomes rigid
                                for i = 1:L
                                    for j = 1:L
                                        A_trial(i,j) = max([A_temp(i+1,j+1),...
                                            A_temp(i-1+1,j+1)*B_temp(i-1+1,j-1+1)*B_temp(i-1+1,j+1),...
                                            A_temp(i+1,j+1+1)*B_temp(i+1,j+1)*B_temp(i-1+1,j+1),...
                                            A_temp(i+1+1,j+1)*B_temp(i+1,j+1)*B_temp(i+1,j-1+1),...
                                            A_temp(i+1,j-1+1)*B_temp(i-1+1,j-1+1)*B_temp(i+1,j-1+1)]);
                                    end
                                end

                                N_rigid_trial(kkk) = sum(sum(A_trial))+sum(sum(B_trial));
                            end

                            switch rule 
                                case 1
                                    [~,kkk_opt] = max(N_rigid_trial); % most efficient 
                                    kkk_opt = kkk_opt(1); 
                                case 2
                                    [~,kkk_opt] = min(N_rigid_trial); % least efficient 
                                    kkk_opt = kkk_opt(1); 
                            end
                        else
                            kkk_opt = 1;
                        end

                        %% update the pattern using the optimal choice

                        chosen = [chosen, ids_batch(kkk_opt)];
                      
                        % exclude the chosen tile from the candidate list  
                        candidates(candidates==ids_batch(kkk_opt)) = [];

                        if ids_batch(kkk_opt) <= L^2
                            A(ids_batch(kkk_opt)) = 1;
                        else
                            B(ids_batch(kkk_opt)-L^2) = 1;
                        end

                        A_temp = zeros(L+2,L+2);
                        B_temp = zeros(L+1,L+1);
                        A_temp(2:end-1,2:end-1) = A;
                        B_temp(2:end-1,2:end-1) = B;

                        % check if a hole becomes rigid
                        for i = 1:(L-1)
                            for j = 1:(L-1)
                                B(i,j) = max([B_temp(i+1,j+1),...
                                    B_temp(i,j+1)*A_temp(i+1,j+1)*A_temp(i+1,j+1+1),...
                                    B_temp(i+1+1,j+1)*A_temp(i+1+1,j+1)*A_temp(i+1+1,j+1+1),...
                                    B_temp(i+1,j-1+1)*A_temp(i+1,j+1)*A_temp(i+1+1,j+1),...
                                    B_temp(i+1,j+1+1)*A_temp(i+1,j+1+1)*A_temp(i+1+1,j+1+1)]);
                            end
                        end

                        % check if a quad becomes rigid
                        for i = 1:L
                            for j = 1:L
                                A(i,j) = max([A_temp(i+1,j+1),...
                                    A_temp(i-1+1,j+1)*B_temp(i-1+1,j-1+1)*B_temp(i-1+1,j+1),...
                                    A_temp(i+1,j+1+1)*B_temp(i+1,j+1)*B_temp(i-1+1,j+1),...
                                    A_temp(i+1+1,j+1)*B_temp(i+1,j+1)*B_temp(i+1,j-1+1),...
                                    A_temp(i+1,j-1+1)*B_temp(i-1+1,j-1+1)*B_temp(i+1,j-1+1)]);
                            end
                        end
                        
                        % total rigid tile and hole count
                        N_rigid = sum(sum(A))+sum(sum(B));
                        N_rigid_all(iii,jjj) = N_rigid;

                        disp(['Rule = ', num2str(rule), ', L = ', num2str(L),...
                            ', k = ', num2str(k),', Sim# = ',num2str(jjj), ...
                            ', Step# = ', num2str(iii), ', #Rigid = ', num2str(N_rigid)]);
                        
                    else
                        % already rigid, can skip the remaining steps
                        N_rigid_all(iii,jjj) = N_max;
                    end
                end

            end

            % save result
            ppsave(L, k, rule, N_rigid_all, n_sim);

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis

cmap = lines(length(L_all));

%% Plot r vs P for each rule and each k

for rule = rule_all
    for k = k_all
        figure;
        hold on;
        for L = L_all

            id = find(L_all == L);

            N_max = L^2+(L-1)^2;
            
            if k == 1
                load(['angle_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['angle_L_',num2str(L),'_k_',num2str(k),'_rule_',num2str(rule),'.mat']);
            end
            
            % percentage of components explicitly rigidified
            r = linspace(0,1,(size(N_rigid_all,1)+1))';
            
            % probability of getting a rigid pattern
            P = sum([zeros(1,n_sim);N_rigid_all]==N_max,2)/n_sim;
            
            plot(r,P,'Color',cmap(id,:),'LineWidth',3);

        end

        legend_cell = cell(1,length(L_all));
        for i = 1:length(L_all)
            legend_cell{i} = ['L = ', num2str(L_all(i))];
        end
        legend(legend_cell);
        title(['Rule = ', num2str(rule),', k = ', num2str(k)]);
        xlabel('r');
        ylabel('P');
        set(gca,'FontSize',16);
        set(gca,'LineWidth',2);
    end
end

%% Plot r vs N/N_max for each rule and each k

for rule = rule_all
    for k = k_all
        figure;
        hold on;
        for L = L_all

            id = find(L_all == L);
            N_max = L^2+(L-1)^2;

            if k == 1
                load(['angle_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['angle_L_',num2str(L),'_k_',num2str(k),'_rule_',num2str(rule),'.mat']);
            end
            
            plot(repmat(linspace(0,1,(size(N_rigid_all,1)+1))',1,n_sim),...
                [zeros(1,n_sim);N_rigid_all]/N_max,'Color',[cmap(id,:),0.1],'LineWidth',1);
        end
        title(['Rule = ', num2str(rule),', k = ', num2str(k)]);
        xlabel('r');
        ylabel('N/N_{max}');
        set(gca,'FontSize',16);
        set(gca,'LineWidth',2);
    end
end

%% Find the critical r for different (L, k)

rc1_all = zeros(length(L_all),length(k_all)); % critical r for rule 1
rc2_all = zeros(length(L_all),length(k_all)); % critical r for rule 2

for rule = rule_all
    for L = L_all
        for k = k_all
            %%
            id = find(L_all == L);
            id2 = find(k_all == k);

            N_max = L^2+(L-1)^2;

            if k == 1
                load(['angle_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['angle_L_',num2str(L),'_k_',num2str(k),'_rule_',num2str(rule),'.mat']);
            end
            
            % probability of getting a rigid pattern
            P = sum([zeros(1,n_sim);N_rigid_all]==N_max,2)/n_sim;

            switch rule
                case 1
                    % critical r: P = 1/2
                    [id_temp] = find(P==0.5);
                    if ~isempty(id_temp)
                        rc1_all(id,id2) = mean(id_temp)/(N_max+1);
                    else
                        rc1_all(id,id2) = interp1(P(find(P<0.5,1,'last'):find(P>0.5,1)),...
                            find(P<0.5,1,'last'):find(P>0.5,1),1/2)/(N_max+1);
                    end

                case 2
                    % critical r: P = (1+plateau)/2
                    if k == 1
                        opt = 1/2;
                    else
                        opt = (1+P(round(length(P)*0.8)))/2;
                    end

                    [id_temp] = find(P==opt);
                    if ~isempty(id_temp)
                        rc2_all(id,id2) = mean(id_temp)/(N_max+1);
                    else
                        rc2_all(id,id2) = interp1(P(find(P<opt,1,'last'):find(P>opt,1)),...
                            find(P<opt,1,'last'):find(P>opt,1),opt)/(N_max+1);
                    end
            end
        end
    end
end

% log-log plot (L = 20) 
L = 20;
id = find(L_all == L);

% Rule 1, lower bound = (4*L-3)/(L^2+(L-1)^2)
figure;
plot(log(k_all),log(rc1_all(id,:)-(4*L-3)/(L^2+(L-1)^2)),'o',...
    'MarkerFaceColor',[201,0,22]/255,'Color',[201,0,22]/255);
lsline
xlabel('log k');
ylabel('log (r_c - r_{min})');
title('Rule 1');
set(gca,'FontSize',16);
set(gca,'LineWidth',2);

% Rule 2, upper bound = 1
figure;
plot(log(k_all),log(1-rc2_all(id,:)),'o','MarkerFaceColor',[201 0 22]/255,...
    'Color',[201,0,22]/255);
xlabel('log k');
ylabel('log (1-r_c)');
title('Rule 2');
set(gca,'FontSize',16);
set(gca,'LineWidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the simulation results
function ppsave(L,k,rule,N_rigid_all,n_sim) 
    if k == 1
        % random, no need to distinguish between the rules
        save(['result_angle/angle_L_',num2str(L),'_k_',num2str(k),'.mat'],...
            'L','N_rigid_all','k','n_sim');
    else
        save(['result_angle/angle_L_',num2str(L),'_k_',num2str(k),'_rule_',...
            num2str(rule),'.mat'],'L','N_rigid_all','k','n_sim');
    end
end
