function d = changeVolatility(d,v)
% function d = changeVolatility(d)
% 
% Change volatility of probability changes
% d = dataset
% v = proportion of volatility change from existing

oldP = d.P;
oldEV = d.EV;
nTrls = size(d,1);

% Smooth over probability changes
oldVolatility = var(diff(oldP));
targetVolatility = v*oldVolatility;

roundVolatility = 3;

currentVolatility = oldVolatility;
switchProb = 0.5;
while true

    % Generate probabilities
    switchProb = switchProb - 0.01;
    if switchProb<=0 || switchProb >=1
        switchProb = 0.5;
    end

    d.P = nan(nTrls,1);
    for trl = 1:nTrls
        if trl==1
            d.P(trl) = datasample([0.1 0.3 0.5 0.7 0.9],1);
        else
            prevProb = d.P(trl-1);
            if datasample([true false],1,'weights',[switchProb 1-switchProb])
                d.P(trl) = datasample(setdiff([0.1 0.3 0.5 0.7 0.9],prevProb),1);
            else
                d.P(trl) = prevProb;
            end
        end
    end

    currentVolatility = var(diff(d.P));

    disp(['Current = ' num2str(round(currentVolatility,roundVolatility)) ', target = ' num2str(round(targetVolatility,1))])
    
    % Recompute EV
    d.EV = d.P.*d.nV_1 + (1-d.P).*d.nV_2;
    h = kstest2(oldEV,d.EV);

    % Do checks
    if round(currentVolatility,roundVolatility) == round(targetVolatility,roundVolatility) && h==0
        break
    end

end

% Plot change
figure
subplot(3,1,1)
plot(oldP); hold on; plot(d.P); xlabel('Trials'); legend({'Before','After'}); ylabel('Probability')
subplot(3,1,2)
histogram(oldEV); hold on; histogram(d.EV); xlabel('EV'); legend({'Before','After'})
subplot(3,1,3)
bar([mean(oldEV>0) mean(oldEV<0); mean(d.EV>0) mean(d.EV<0)]')
legend({'Before','After'}); title('Correct choices'); set(gca,'xticklabels',{'Approach','Avoid'})

end