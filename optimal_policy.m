%% Compute optimal policy

clear all
close all
clc

% Action/map parameters
actions = 1:4;
action_labels = {'UP','DOWN','LEFT','RIGHT'};
response_text = {'^','v','<','>'};

map = [0 0 7 0 4
       0 0 6 0 3
       9 8 5 1 2];

risk_state = [5 6];
predator_state = 7;
reward_state = 9;
neutral_state = 4;

start_state = 1;
start_health = [.4 .5 .6];

% Task parameters
risk_prob = [0 .15 .3];
risk_names = {'ASLEEP','HALF-AWAKE','ALERT'};
reward_mag = .2;
attack_mag = .5; % certain attack value OR set to Inf so that you die on attack
starve_mag = .1;

disp(['R = ' num2str(reward_mag*100) '%'])
disp(['L = ' num2str(attack_mag*100) '%'])

% small_reward_mag = .1;
% small_reward_prob = 0;

nDays = 5; % how many days per trial (in both practice and real task)
     
%% Transition matrix between health 'states'

health_idx = round(0:.1:1,1);
actions = {'safe','risky'}; % let's assume that these are the only 2 choices you're really making

TM = zeros(length(risk_prob),length(actions),length(health_idx),length(health_idx));
RM = zeros(length(risk_prob),length(actions),length(health_idx),length(health_idx));
figure
cc = 0;
for r = 1:length(risk_prob)
    for a = 1:length(actions)
        for h1 = 1:length(health_idx) % health BEFORE
            
            if health_idx(h1) == 0 % if dead
                TM(r,a,h1,1) = 1; % a transition from a dead state to another dead state is given a 1
                
            else % if not dead
                
                if a == 1 % no reward, just starve
                    
                    h2 = health_idx(h1) - starve_mag;
                    if h2 < 0
                        h2 = 0;
                    end
                    h2 = round(h2,1);
                    TM(r,a,h1,health_idx == h2) = 1; % with certainty
                    
                elseif a == 2 % conflict path
                    
                    % get past predator
                    h2 = health_idx(h1) + reward_mag;
                    if h2 > 1
                        h2 = 1;
                    end
                    h2 = round(h2,1);
                    TM(r,a,h1,health_idx == h2) = 1-risk_prob(r);
                    
                    % get attacked
                    h2 = health_idx(h1) - attack_mag;
                    if h2 < 0
                        h2 = 0;
                    end
                    TM(r,a,h1,health_idx == h2) = risk_prob(r);
                    
                end
            end
        end
        tmp = squeeze(TM(r,a,:,:));
        TM(r,a,:,:) = tmp./repmat(nansum(tmp,2),1,size(tmp,1)); % normalise probabilities to 1
        
        RM(r,a,:,1) = -1; % all transitions to zero energy state (column 1) give penalty of -1
        RM(r,a,:,end) = 1; % all transitions to maximum energy state (column 1) give reward of 1
        
        cc = cc + 1;
        subplot(length(risk_prob),length(actions),cc)
        imagesc(squeeze(TM(r,a,:,:)))
        set(gca,'XTick',1:length(health_idx))
        set(gca,'YTick',1:length(health_idx))
        set(gca,'XTickLabels',{'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'})
        set(gca,'YTickLabels',{'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'})
        title([risk_names{r} ', chose ' actions{a}])
        ylabel('start health')
        xlabel('end health')
        colormap('hot')
        hold on
        plot([0 length(health_idx)],[0 length(health_idx)],'g')
        colorbar
        
    end
end
    

% Reward at the end of the week, depending on final health
terminal_reward = ones(1,length(health_idx));
terminal_reward(end) = 2;
terminal_reward(1) = 0;

%% Optimal Policy (Backward Induction)

N = nDays;
discount = 1;

opt_policy = {};

% figure options
cmap = [0 0 0 % blacked out
        0 0 1    % go hungry (blue)
        1 0 0    % conflict path (red)
        .5 0 .5];    % either one (purple)

choice_summary = nan(length(risk_prob),3);
for r = 1:length(risk_prob)

    P = permute(squeeze(TM(r,:,:,:)),[2 3 1]);
    R = permute(squeeze(RM(r,:,:,:)),[2 3 1]);
    
    PR = []; % compute the probabilistic reward for being in each state (row) depending on each choice (column: starve, risky)
    for a = 1:size(P,3) % for the number of actions...
        PR(:,a) = sum(P(:,:,a).*R(:,:,a),2);
    end
    
    S = size(P,1); % no. of states
    V = zeros(S,N+1); % value function
    V(:,N+1) = terminal_reward;
    
    policy = []; % optimal choice at each health (row) given the day (col)
    for n = 0:N-1
        
        Vprev = V(:,N-n+1);
        
        Q = []; % row = states (health), col = action (safe vs. risk path)
        for a = 1:size(P,3)
            Q(:,a) = PR(:,a) + discount * P(:,:,a) * Vprev;
        end
        [W, X] = max(Q,[],2);
        X(Q(:,1) == Q(:,2)) = 3; % 3 = indifferent
        
        V(:,N-n) = W; % the maximum Q values
        policy(:,N-n) = X; % the actions that had maximum Q values

    end
    
    opt_policy{r} = policy;

    figure
    fig_pol = policy;
    fig_pol(round(setdiff(health_idx,start_health)*10)+1,1) = 0; % black out non-starting health states
    for d = 2:size(policy,2)
        prev_health = health_idx(fig_pol(:,d-1)>0);
        fig_pol(round(setdiff(health_idx,round([0 prev_health+reward_mag, prev_health-starve_mag],1))*10)+1,d) = 0; % black out impossible subsequent states (given reward/loss mags)
    end
    
    choice_prob(r,1) = sum(fig_pol(:) == 1) / sum(fig_pol(:) ~= 0);
    choice_prob(r,2) = sum(fig_pol(:) == 2) / sum(fig_pol(:) ~= 0);
    choice_prob(r,3) = sum(fig_pol(:) == 3) / sum(fig_pol(:) ~= 0);
    
    imagesc(fig_pol(2:end,:))
    colormap(cmap)
    title([num2str(risk_prob(r)) '% predator probability'])
    set(gca,'XTick',1:N)
    set(gca,'YTick',1:length(health_idx))
    set(gca,'YTickLabels',{'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'})
    xlabel('Day')
    ylabel('Health')
    
end

figure
bar(choice_prob);
set(gca,'XTickLabels',risk_prob)
xlabel('P(Predator)')
ylabel('% Choices across days')
legend({'Go Hungry','Conflict Path','Whichever'})
title('Choice Distribution across Health/Days')