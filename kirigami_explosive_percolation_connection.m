% Simulation for explosive rigidity percolation using the connection-based method
% 
% - Add a connection between two tiles each step based on some selection rule
% - Check the total DOF of the entire pattern
%
% Reference:
% G. P. T. Choi, L. Liu, L. Mahadevan, "Explosive rigidity perolcation in
% kirigami", preprint, 2022.
% 
% Copyright (c) 2022,  Gary P. T. Choi, L. Liu, L. Mahadevan


mkdir('result_connection');
addpath('result_connection');

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
        
        n_maxlink = 4*L*(L-1);

        mat = zeros(12*L^2,3);

        % Edge length constraints
        % 4 quad boundary constraints, and 1 no shear constraints (Direction fixed
        % for now, from bottom left to top right)
        for i = 1:L^2
            mat(i*12-11,:)=[i*5-4,i*8-7,-1];
            mat(i*12-10,:)=[i*5-4,i*8-5,1];
            mat(i*12-9,:)=[i*5-3,i*8-7,1];
            mat(i*12-8,:)=[i*5-3,i*8-3,-1];
            mat(i*12-7,:)=[i*5-3,i*8-6,1];
            mat(i*12-6,:)=[i*5-3,i*8-2,-1];
            mat(i*12-5,:)=[i*5-2,i*8-6,1];
            mat(i*12-4,:)=[i*5-2,i*8-0,-1];
            mat(i*12-3,:)=[i*5-1,i*8-4,1];
            mat(i*12-2,:)=[i*5-1,i*8-2,-1];
            mat(i*12-1,:)=[i*5-0,i*8-3,1];
            mat(i*12-0,:)=[i*5-0,i*8-1,-1];
        end

        mat_ori = mat;

        % Construct the list of all possible connections
        linkpairs = [];

        % Boundary connections
        for i = 1:(L-1)
            linkpairs(end+1,:)=[4*i-2,4*(i+1)-3];
        end
        for i = (L^2-L+1):(L^2-1)
            linkpairs(end+1,:)=[4*i-1,4*(i+1)];
        end
        for i = 1:L:(L^2-L)
            linkpairs(end+1,:)=[4*i,4*(i+L)-3];
        end
        for i = L:L:(L^2-L)
            linkpairs(end+1,:)=[4*i-1,4*(i+L)-2];
        end

        % Inner connections
        % horizontal
        for jj = 1:(L-1)
            for i = (jj*L-(L-1)):(jj*L-1)
                linkpairs(end+1,:)=[4*i-1,4*(i+1)];
            end
        end
        for jj = 1:(L-1)
            for i = (jj*L+1):(jj*L+L-1)
                linkpairs(end+1,:)=[4*i-2,4*(i+1)-3];
            end
        end
        % vertical
        for jj = 1:(L-1)
            for i = jj:L:(L^2-L)
                linkpairs(end+1,:)=[4*i-1,4*(i+L)-2];
            end
        end
        for jj=1:L-1
            for i = (jj+1):L:(L^2-L)
                linkpairs(end+1,:)=[4*i,4*(i+L)-3];
            end
        end

        for k = k_all
            dof_all = zeros(n_maxlink,n_sim);
            parfor jjj = 1:n_sim
                %% start a simulation

                mat = mat_ori;
                candidates = 1:n_maxlink;
                chosen = [];
                
                dof_old = 3*L^2;
                
                for iii = 1:n_maxlink
                    
                    if dof_old ~= 3
                        % k choices
                        if length(candidates) ~= 1
                            ids_batch = randsample(candidates,min(length(candidates),k));
                        else
                            ids_batch = candidates;
                        end
                        if k > 1
                            %% consider the outcome for each choice
                            dof_temp = zeros(min(length(candidates),k),1);
                            for kkk = 1:min(length(candidates),k)
                                % Add the connection
                                newmat = [mat; ... 
                                          max(mat(:,1))+1,linkpairs(ids_batch(kkk),1)*2-1,1;...
                                          max(mat(:,1))+1,linkpairs(ids_batch(kkk),2)*2-1,-1; ...
                                          max(mat(:,1))+2,linkpairs(ids_batch(kkk),1)*2,1; ...
                                          max(mat(:,1))+2,linkpairs(ids_batch(kkk),2)*2,-1];

                                % Calculate the DoF (temporary)
                                rgd_Matrix = sparse(newmat(:,1),newmat(:,2),newmat(:,3),...
                                    size(mat,1)+2, 8*L^2);
                                [r] = calc_rank(rgd_Matrix);
                                dof_temp(kkk) = 8*L^2-r;

                            end

                            switch rule 
                                case 1
                                    [~,kkk_opt] = min(dof_temp); % most efficient
                                    kkk_opt = kkk_opt(1);
                                case 2
                                    [~,kkk_opt] = max(dof_temp); % least efficient
                                    kkk_opt = kkk_opt(1);
                            end
                        else
                            kkk_opt = 1;
                        end

                        %% update the pattern using the optimal choice
                        
                        chosen = [chosen, ids_batch(kkk_opt)];

                        newmat = [mat; ... 
                                  max(mat(:,1))+1,linkpairs(ids_batch(kkk_opt),1)*2-1,1;...
                                  max(mat(:,1))+1,linkpairs(ids_batch(kkk_opt),2)*2-1,-1; ...
                                  max(mat(:,1))+2,linkpairs(ids_batch(kkk_opt),1)*2,1; ...
                                  max(mat(:,1))+2,linkpairs(ids_batch(kkk_opt),2)*2,-1];

                        % Calculate the DoF
                        rgd_Matrix = sparse(newmat(:,1),newmat(:,2), ...
                            newmat(:,3), size(mat,1)+2, 8*L^2);
                        [r]=calc_rank(rgd_Matrix);
                        dof = 8*L^2-r;
                        dof_all(iii,jjj) = dof;

                        dof_old = dof;
                        mat = newmat;
                        
                        % exclude the chosen connection from the candidate list
                        candidates = setdiff(candidates,ids_batch(kkk_opt));

                        disp(['Rule = ', num2str(rule), ', L = ', num2str(L),...
                            ', k = ', num2str(k),', Sim# = ',num2str(jjj), ...
                            ', Step# = ', num2str(iii), ', #DOF = ', num2str(dof)]);
                        
                    else
                        % already rigid, can skip the remaining steps
                        dof_all(iii,jjj) = 3;
                    end
                end

            end

            % save result
            ppsave(L, k, rule, dof_all, n_sim);

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
        
            if k == 1
                load(['connection_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['connection_L_',num2str(L),'_k_',num2str(k),...
                    '_rule_', num2str(rule),'.mat']);
            end
            
            % percentage of connections added
            r = linspace(0,1,(size(dof_all,1)+1))';
            
            % probability of getting a rigid pattern
            P = sum([3*L^2*ones(1,n_sim);dof_all]==3,2)/n_sim;

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


%% Plot r vs (DOF-3)/(3L^2-3) for each rule and each k

for rule = rule_all
    for k = k_all
        figure;
        hold on;
        for L = L_all
            
            id = find(L_all == L);
    
            if k == 1
                load(['connection_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['connection_L_',num2str(L),'_k_',num2str(k),...
                    '_rule_', num2str(rule),'.mat']);
            end
            
            plot(repmat(linspace(0,1,(size(dof_all,1)+1))',1,n_sim),...
                [3*L^2*ones(1,n_sim)-3;dof_all-3]/(3*L^2-3),'Color',...
                [cmap(id,:),0.1], 'LineWidth',1);
        end
        title(['Rule = ', num2str(rule),', k = ', num2str(k)]);
        xlabel('r');
        ylabel('(d-3)/(3L^2-3)');
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

            n_maxlink = 4*L*(L-1);

            if k == 1
                load(['connection_L_',num2str(L),'_k_',num2str(k),'.mat']);
            else
                load(['connection_L_',num2str(L),'_k_',num2str(k),...
                    '_rule_', num2str(rule),'.mat']);
            end

            % probability of getting a rigid pattern
            P = sum([3*L^2*ones(1,n_sim);dof_all]==3,2)/n_sim;

            switch rule
                case 1
                    % critical r: P = 1/2
                    [id_temp] = find(P==0.5);
                    if ~isempty(id_temp)
                        rc1_all(id,id2) = mean(id_temp)/(n_maxlink+1);
                    else
                        rc1_all(id,id2) = interp1(P(find(P<0.5,1,'last'):find(P>0.5,1)),...
                            find(P<0.5,1,'last'):find(P>0.5,1),1/2)/(n_maxlink+1);
                    end
                case 2
                    % critical r: P = 1/2
                    [id_temp] = find(P==0.5);
                    if ~isempty(id_temp)
                        rc2_all(id,id2) = mean(id_temp)/(n_maxlink+1);
                    else
                        rc2_all(id,id2) = interp1(P(find(P<0.5,1,'last'):find(P>0.5,1)),...
                            find(P<0.5,1,'last'):find(P>0.5,1),1/2)/(n_maxlink+1);
                    end
            end
        end
    end
end

% log-log plot (L = 20) 
L = 20;
id = find(L_all == L);

% Rule 1, lower bound = ceil((3*L^2-3)/2)/(4*L*(L-1))
figure;
plot(log(k_all),log(rc1_all(4,:)-ceil((3*L^2-3)/2)/(4*L*(L-1))),'o',...
    'MarkerFaceColor',[201,0,22]/255,'Color',[201,0,22]/255);
lsline
xlabel('log k');
ylabel('log (r_c - r_{min})');
title('Rule 1');
set(gca,'FontSize',16);
set(gca,'LineWidth',2);

% Rule 2, upper bound = 1
figure;
plot(log(k_all), log(1-rc2_all(id,:)),'o','MarkerFaceColor',[201,0,22]/255,...
    'Color',[201,0,22]/255);
xlabel('log k');
ylabel('log (1-r_c)');
title('Rule 2');
set(gca,'FontSize',16);
set(gca,'LineWidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the simulation results
function ppsave(L,k,rule,dof_all,n_sim) 
    if k == 1
        % random, no need to distinguish between the rules
        save(['result_connection/connection_L_',num2str(L),'_k_',num2str(k),...
            '.mat'],'L','dof_all','k','n_sim');
    else 
        save(['result_connection/connection_L_',num2str(L),'_k_',num2str(k),...
            '_rule_',num2str(rule),'.mat'],'L','dof_all','k','n_sim');
    end
end