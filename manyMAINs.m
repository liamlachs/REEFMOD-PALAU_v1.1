%__________________________________________________________________________
%
% Local scheduling script
%
% This script is used to submit different REEFMOD instances locally, with specific parameters
%
% Liam Lachs, liamlachs@gmail.com, 12/2023
%__________________________________________________________________________


% for gcm = 12:16 % 11 %[11 13 14] % 1:16 %[2 5 6 7 8 10 13 14 16]
%     for ssp = [ 1 2 4] 
%         for h2 = 0:0.1:1 
%             MAIN_REEFMOD_PAL(h2,ssp,gcm,0,0,0)
%         end
%     end
% end


% parfor gcm = [1:16] 
%      MAIN_REEFMOD_PAL(0.3, 1  ,gcm ,0,0,0)
% end

parfor gcm = [1:16] % 
%     for ssp = [ 1 2 4 ] 
        for ET_multiplier = [-0.5:.1:-0.1, 0.1:.1:-0.5]% [-1, -0.8 -0.6, 0.6, 0.8, 1] % [-0.4, -0.2, 0.2, 0.4]
            m = ET_multiplier;
            MAIN_REEFMOD_PAL(0.3, 1  ,gcm ,m,m,-m)
        end
%     end
end