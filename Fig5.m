%Script to simulate data from Fig 5 of Verharen et al., Psychopharmacology
%2020.
%Written by Jeroen P.H. Verharen, 2019.
%Just click 'Run', and it should work! (Tested in Matlab 2018b on Mac OS)

%% Main script - loop through different {a+, a-, pi} conditions, and arrange in 3-dimensional grid

xrange = 0:0.02:1; %x-axis range (reward learning)
yrange = 0:0.02:1; %y-axis range (punishment learning)
zrange = [-0.25:0.25:0.75]; %z-axis range (stickiness)

%pre-allocate vectors to save convential performance measures from
%simulated data
matRev = nan(length(xrange), length(yrange), length(zrange));
matRew = nan(length(xrange), length(yrange), length(zrange));
matWin = nan(length(xrange), length(yrange), length(zrange));
matLose = nan(length(xrange), length(yrange), length(zrange));

for x = 1:length(xrange) %loop through reward learning
    
    fprintf('Progress: %.0f of %.0f...\n', x, length(xrange)) %print progress
    
    for y = 1:length(yrange) %loop through punishment learning
        for z = 1:length(zrange) %loop through stickiness
            
            vecOutput = revlearn_simulate(xrange(x), yrange(y), zrange(z), 50); %run 50 simulation
            
            %collect data and arrange in 3D grid
            matRev(x,y,z) = mean(vecOutput.reversals);
            matRew(x,y,z) = mean(vecOutput.rewarded);
            matWin(x,y,z) = mean(vecOutput.winstay);
            matLose(x,y,z) = mean(vecOutput.loseshift);
            
        end
    end
end

%% plot
figure
colormap(parula)

for cols = 1:4 %loop through columns (types of performance measure)

    if cols == 1
        toplot = imgaussfilt(matRev, .75);
    elseif cols == 2
        toplot = imgaussfilt(matRew, .75);
    elseif cols == 3
        toplot = imgaussfilt(matWin, .75);
    elseif cols == 4
        toplot = imgaussfilt(matLose, .75);
    end

for k = 1:5 %loop through 5 steps of z dimension (stickiness values)
    subplot(5,4,(k-1)*4+cols)
    imagesc(toplot(:,:,k), [min(toplot(:)) max(toplot(:))])
    hold on
    contour(toplot(:,:,k),  [prctile(toplot(:), 90) prctile(toplot(:), 90)], 'LineColor', 'k')
    
    title(sprintf('pi = %.2f', zrange(k)))
    xlabel('a+')
    ylabel('a-')
    set(gca,'XTick', [1 length(xrange)], 'XTickLabel', [0 1])
    set(gca,'YTick', [1 length(xrange)], 'YTickLabel', [0 1])
    colorbar
end

end











%% function containing simulator

function vecOutput = revlearn_simulate(dblLearningRateUpval, dblLearningRateDeval, pi, iterations)
%% set params

%% pre allocate vectors for conventional
vecOutput.reversals = nan(iterations,1);
vecOutput.rewarded = nan(iterations,1);
vecOutput.winstay = nan(iterations,1);
vecOutput.loseshift = nan(iterations,1);


for iter = 1:iterations
    
    %% pre allocate vectors for modeling
    vecChoice = nan(200,1);
    vecActive = nan(200,1);
    vecWinLose = nan(200,1);
    
    
    %% run simulation
    
    
    boolLeft = 0;
    boolRight = 0;
    valLeft = 0;
    valRight = 0;
    dblBeta = 1.686; %value of beta fixed at grand average of Fig. 3
    active = 1;
    activesum = 0;
    totalreversals = 0;
    
    lowOdds = [1 0 0 0 0]; %odds of winning low probability hole
    highOdds = [0 1 1 1 1]; %odds of winning high probability hole
    
    for i = 1 : 200 %loop through trials
        
        pRight = exp(dblBeta*valRight + pi*boolRight)/( exp(dblBeta * valLeft + pi*boolLeft) + exp(dblBeta * valRight + pi*boolRight)); %softmax
        
        %make choice
        if pRight > rand %right choice
            choice = 1;
        else
            choice = 0;
        end
        
        %task is probabilistic, so determine if reward is obtained in this
        %trial
        if choice == active %if high probabily hole is chosen
            activesum = activesum + 1;
            winlose = randsample(highOdds,1);
        else % if low probability hole is chosen
            winlose = randsample(lowOdds,1);
            activesum = 0;
        end
        
        %update value
        if choice == 1 &&  winlose == 1 %right chosen, right won
            valRight = valRight + dblLearningRateUpval*(1 - valRight); %Rescorla-Wagner
        elseif choice == 1 &&  winlose == 0 %right chosen, right lost
            valRight = valRight + dblLearningRateDeval*(0 - valRight); %Rescorla-Wagner
        elseif choice == 0 &&  winlose == 1 %left chosen, left won
            valLeft = valLeft + dblLearningRateUpval*(1-valLeft); %Rescorla-Wagner
        elseif choice == 0 &&  winlose == 0 %left chosen, left lost
            valLeft = valLeft + dblLearningRateDeval*(0-valLeft); %Rescorla-Wagner
        end
        
        %set booleans for stickiness
        if choice == 1
            boolRight = 1;
            boolLeft = 0;
        elseif choice == 0
            boolLeft = 1;
            boolRight = 0;
        end
        
        %do reversal after 8 choices for high=prob
        if activesum == 8
            
            if active == 1
                active = 0;
            elseif active == 0
                active = 1;
            end
            
            activesum = 0;
            totalreversals = totalreversals + 1;
            
        end
        
        %save trial-to-trial data
        vecChoice(i) = choice;
        vecActive(i) = active;
        vecWinLose(i) = winlose;
        
        
    end
    
    %% analyse behavioral data
    
    totalwin = 0;
    totalwinstay = 0;
    
    totallose = 0;
    totalloseshift = 0;
    
    for i = 1: 199
        
        if vecWinLose(i) == 1 %determine behavior in win trials
            
            totalwin = totalwin + 1;
            
            if vecChoice(i) == vecChoice(i+1)
                totalwinstay = totalwinstay + 1;
            end
            
        elseif vecWinLose(i) == 0 %determine behavior in lose trials
            
            totallose = totallose + 1;
            
            if vecChoice(i) ~= vecChoice(i+1)
                totalloseshift = totalloseshift + 1;
            end
            
        end
    end
    
    vecOutput.reversals(iter) = 100*totalreversals/200; %compute number of reversals per 100 trials and save data
    vecOutput.rewarded(iter) = sum(vecWinLose)/200; %compute fraction of rewarded trials and save data
    vecOutput.winstay(iter) = totalwinstay/totalwin; %compute win-stay and save data
    vecOutput.loseshift(iter) = totalloseshift/totallose; %compute lose-shift and save data
    
end
end