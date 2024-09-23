%__________________________________________________________________________
%
% Local scheduling script 
%
% This script is used to submit different REEFMOD instances locally, with specific parameters
%
% Liam Lachs, liamlachs@gmail.com, 12/2023
%__________________________________________________________________________

% IPCC best estimate ECS models are:
% [2 3 11 12 13 15]

% 'Hot' models with ECS > 4C are:
% [1 4:10 14 16]

% parfor gcm = 1:16
%         for ssp = [1 2 4]
%             for h2 = 0:0.1:1
%                 MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
%             end
%         end
% end

% i can now utilise 24 cores simulataneously, so i need to call two vectors
% job = batch('manyMAINs','Pool',23)

% that maximise that output. I want to do all best estimate models first
% h2vals = repmat(0:0.1:1,1,3);
h2vals = repmat([0 0.3 1],1,3);
sspvals = repelem([1 2 4],1,length(h2vals)/3);

gcmIPCCbest = [2 3 11 12 13 15];
parfor i = 1:length(h2vals)
    h2 = h2vals(i);
    ssp = sspvals(i);
    for gcm = gcmIPCCbest
        MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
    end
end

gcmHotModels = [1 4:10 14 16];
parfor i = 1:length(h2vals)
    h2 = h2vals(i);
    ssp = sspvals(i);
    for gcm = gcmHotModels
        MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
    end
end

% do only h2=0.3, for the 20 simulations run or other runs that only do h2=3
h2=0.3;
gcmIPCCbest = [2 3 11 12 13 15];
sspvals = [repmat(1,1,length(gcmIPCCbest)) repmat(2,1,length(gcmIPCCbest)) repmat(4,1,length(gcmIPCCbest))];
gcmvals = [repmat(gcmIPCCbest,1,3)];
parfor i = 1:length(sspvals)
    gcm = gcmvals(i);
    ssp = sspvals(i);
    MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
end

gcmHotModels = [1 4:10 14 16];
sspvals = [repmat(1,1,length(gcmHotModels)) repmat(2,1,length(gcmHotModels)) repmat(4,1,length(gcmHotModels))];
gcmvals = [repmat(gcmHotModels,1,3)];
parfor i = 1:length(sspvals)
    gcm = gcmvals(i);
    ssp = sspvals(i);
    MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
end
% parfor gcm = 1:16 % 1:16 % 11 %[11 13 14] % 1:16
% %     if sum([1 4:10 14 16]==gcm)>0
%         for ssp = [ 1 2 4]
%             for h2 = [0.3 1.0] % 0:0.1:1
%                 MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
%             end
%         end
% %     end
% end

% parfor gcm = [1:16] 
%      MAIN_REEFMOD_PAL(0.3, 1  ,gcm ,0,0,0)
% end

% parfor gcm = [1:16] % 
% %     for ssp = [ 1 2 4 ] 
%         for ET_multiplier = [-0.5:.1:-0.1, 0.1:.1:-0.5]% [-1, -0.8 -0.6, 0.6, 0.8, 1] % [-0.4, -0.2, 0.2, 0.4]
%             m = ET_multiplier;
%             MAIN_REEFMOD_PAL(0.3, 1  ,gcm ,m,m,-m)
%         end
% %     end
% end