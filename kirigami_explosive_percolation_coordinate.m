% Simulation for explosive rigidity percolation using the coordinate-based method
% 
% - Rigidify a tile by fixing its vertex coordinates each step 
%   based on some selection rule
% - Check and update the rigidity of other tiles/holes based on the
%   rigidity propagation rule
% - Count the total number of rigid tiles and holes
%
% Reference:
% G. P. T. Choi, L. Liu, L. Mahadevan, "Explosive rigidity perolcation in
% kirigami", preprint, 2022.
% 
% Copyright (c) 2022,  Gary P. T. Choi, L. Liu, L. Mahadevan

addpath(genpath('result_coordinates'));

%% Setup

% kirigami pattern size
L_all = [5,10,15,20];

% number of choices in each step
k_all = [1,2,5,10,25];

% selection rule (1: most efficient, 2: least efficient)
rule_all = [1,2]; 

% number of simulations for each setup
n_sim = 200;

%% Run the simulations 
% refer to the python code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis

cmap = lines(length(L_all));

%% Plot r vs P for each rule and each k

for rule = rule_all
    if rule == 1
        M = readmatrix('pcts_rigid_grouped_maximize_rigidity.csv');
    else
        M = readmatrix('pcts_rigid_grouped_minimize_rigidity.csv');
    end
    
    for k = k_all
        figure;
        hold on;
        for L = L_all
            id = find(L_all == L);
            idm = find(M(:,2) == L & M(:,5) == k);
            
            % percentage of tiles explicitly rigidified
            r = M(idm,3);
            
            % probability of getting a rigid pattern
            P = M(idm,9);
            
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
    if rule == 1
        M = readmatrix('max_rigidity/ind_sim_results.csv');
    else
        M = readmatrix('min_rigidity/ind_sim_results.csv');
    end
    for k = k_all
        figure;
        hold on;
        for L = L_all
            id = find(L_all == L);
            idm = find(M(:,1) == L & M(:,5) == k);
            N_max = L^2+(L-1)^2;
            
            % percentage of tiles explicitly rigidified
            r = reshape(M(idm,2),L^2,n_sim);
            
            % total rigid tile and hole count
            N = reshape(M(idm,8)+M(idm,9),L^2,n_sim);
            
            plot(r,N/N_max,'Color',[cmap(id,:),0.1],'LineWidth',1);
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
    if rule == 1
        M = readmatrix('pcts_rigid_grouped_maximize_rigidity.csv');
    else
        M = readmatrix('pcts_rigid_grouped_minimize_rigidity.csv');
    end
    
    for L = L_all
        for k = k_all
            %%
            id = find(L_all == L);
            id2 = find(k_all == k);

            idm = find(M(:,2) == L & M(:,5) == k);
            
            % percentage of tiles explicitly rigidified
            r = M(idm,3);
            
            % probability of getting a rigid pattern
            P = M(idm,9);
            
            N_max = L^2;

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

% Rule 1, lower bound = 2L (for even L) or 2L-1 (for odd L)
figure;
plot(log(k_all),log(rc1_all(id,:)-(2*L-mod(L,2))/L^2),'o',...
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
ylim([-6 -1.5])
set(gca,'FontSize',16);
set(gca,'LineWidth',2);

