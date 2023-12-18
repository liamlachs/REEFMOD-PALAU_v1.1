% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
%
% Last modified: 07/2020 (track coral loss per species)
%
% Modified (08/2018) to integrate genetic adaptation to warming
% -------------------------------------------------------------------------


function [coral, genes, algal, total_coral_loss, total_mortality] = f_bleaching_Liam(coral, genes, algal, bleaching_whole_mortality, DHW, ...
    CORAL, doing_3D, nb_coral_types, doing_clades, doing_DHWbleaching, doing_DHWadaptation_covZw, doing_coral_age, doing_genetics, bleaching_whole_offset, bleaching_partial_offset,Topt_baseline, Topt2index)

% NEED as inputs
% - matrix of DHW for the given reef for each MMM (Topt) at time step
% Extract data from the structures (need to be filled again when leaving)
algal_cm2 = [algal.cover_cm2] ;
[coral_cm2, surface_cm2, volume_cm3, coral_age, clade, heat_tolerance, colony_ID, species_ID] = f_struct_deploy (coral);

% Locate brooders and spawners
% id_brooders = zeros(size(coral_cm2)) ; % initialize brooders id
% id_spawners = zeros(size(coral_cm2)) ; % initialize spawners id

%%%% This is new stuff (August 2013) for implementing species-specific bleaching mortalities
id0 = zeros(size(coral_cm2)) ;

sensitivity_bleaching = id0 ;
extent_bleaching = id0 ;
proba_switching = id0 ;
% expand colony fecundity values for all species
% if doing_DHWadaptation_covZw == 1
%     fecund_min_size = id0;
%     fecund_a = id0;
%     fecund_b = id0;
% end

id1 = id0+1;
id1(coral_cm2 <= 0) = 0 ; % Remove the already dead ones

col_start = 1;
col_stop = 0;

%% Added thermal resistance from genetics (August 2018)
if doing_genetics == 1
    
    MORT = id0 ;
    
    for s = 1:nb_coral_types
        
        if species_ID(s)>0
            
            col_stop = col_stop + species_ID(s) ;
            
            id1_tmp = id1(:,col_start:col_stop);
            s_mort = zeros(size(id1_tmp)) ;
            colony_ID_tmp = colony_ID(:,col_start:col_stop);
            
            % First update the list of QTLs to only keep the survivors from previous mortality
            list_old = genes(s).list_coral_ID;
            list_new = colony_ID_tmp(id1_tmp==1);
            
            check = ismember(list_old,list_new) ;
            
            genes(s).QTLs(check==0,:,:)=[];
            genes(s).list_coral_ID(check==0,:)=[];
            genes(s).phenotypes(check==0,:)=[];
            
            % Phenotype is Topt so need to find out the corresponding DHW
            % Then convert into whole colony mortality using Hughes relationship
            % Output must be a matrix of proba of whole colony-mortality
            All_Topts = round(10*(Topt_baseline + genes(s).phenotypes))/10 ;
            
            if isempty(All_Topts)==0
                
%                 max(All_Topts)               
                s_mort(id1_tmp==1) = bleaching_whole_mortality(10*All_Topts - Topt2index) ;
                
            end
            
            sensitivity_bleaching(:,col_start:col_stop) = CORAL.sensitivity_bleaching(s);
            extent_bleaching(:,col_start:col_stop) = CORAL.bleaching_partial_extent(s);
            proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
            MORT(:,col_start:col_stop) = s_mort ;
            
            col_start = col_start + species_ID(s) ;
            
        end
    end
    
elseif doing_DHWbleaching == 1
    
%     MORT = bleaching_whole_mortality.*id1 ;
    
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s);
        % bleaching sensitivity is scaled to corymbose Acropora (value of
        % 1.4 compared to others) as our new bleaching, mortality
        % probability is for a corymbose Acropora.
        sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s) / 1.4;
        extent_bleaching(:,col_start:col_stop)= CORAL.bleaching_partial_extent(s);
        proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
%         if doing_DHWadaptation_covZw == 1
%             fecund_min_size(:,col_start:col_stop) = CORAL.fecund_min_size(s);
%             fecund_a(:,col_start:col_stop) = CORAL.fecund_a(s);
%             fecund_b(:,col_start:col_stop) = CORAL.fecund_b(s);
%         end
        col_start = col_start + species_ID(s) ;
        
    end

else
    
    MORT = bleaching_whole_mortality.*id1 ;
    
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s);
        extent_bleaching(:,col_start:col_stop)= CORAL.bleaching_partial_extent(s);
        proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
        
        col_start = col_start + species_ID(s) ;
        
    end
    
end

if doing_clades == 1
    % Clade-induced tolerance to thermal stress
    sensitivity_bleaching(clade==2) = sensitivity_bleaching(clade==2) * CORAL.bleaching_tolerance_clade ;
end

%________________________________
%
% Whole-colony mortality
%________________________________



% Generate random mortalities for adol + adults
id1 = ones(size(coral_cm2)) ; % assigns 1 to every colony
% id1(coral_cm2 < CORAL.adol_size) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)
id1(coral_cm2 < CORAL.size_threshold_wcm) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)

rand_mort1 = rand(size(id1)) ; % Generates random probability for adol + adults

if doing_DHWbleaching == 1
    
    % make a per colony bleaching mortality probability based on DHW and heat tolerance
    % logit(bleaching mortality) = B0 + B1 * DHW, where B0 is the heat tolerance of a colony
    %
    %                           exp(B0 + B1 * DHW)
    % bleaching mortality  = ------------------------
    %                         1 + exp(B0 + B1 * DHW)
    %
    LinearEq = exp(heat_tolerance + CORAL.DHWbleaching_mortality_Slope * DHW);
    bleaching_mortality_prob =  (LinearEq ./ (1+LinearEq));
    bleaching_mortality_prob(coral_cm2 <= 0) = 0; % remove already dead ones
    
%     scatter(reshape(heat_tolerance,prod(size(heat_tolerance)),1), ...
%         reshape(bleaching_mortality_prob,prod(size(bleaching_mortality_prob)),1))

    % offset bleaching probabilities for other taxa and based on depth
%     prob_whole_mortality = sensitivity_bleaching * CORAL.bleaching_depth .* bleaching_mortality_prob; 
    % offset bleaching probabilities based on depth
    prob_whole_mortality = CORAL.bleaching_depth .* bleaching_mortality_prob; 

%         prob_whole_mortality = 1-(1-prob_whole_mortality).^2;
%     load('prob_whole_mortality.mat', 'prob_whole_mortality');

else
    
    prob_initial_mortality = sensitivity_bleaching*CORAL.bleaching_depth.*MORT;
    prob_initial_mortality(prob_initial_mortality>1)=1; % cap to 1 before extrapolating to 6 month
    prob_whole_mortality = 1-(1-prob_initial_mortality).^bleaching_whole_offset;
    
    % prob_whole_mortality = id1.*(1-(1-sensitivity_bleaching.*MORT).^bleaching_whole_offset);
    % prob_whole_mortality = id1.*sensitivity_bleaching .* MORT * bleaching_whole_offset;
end
unique(prob_whole_mortality');

% % Compute the response to selection if doing the Price Equation approach to
% % modelling adaptation. Only of gravid colonies.
% % Divide by the generation time (g, not implemented yet)
% % R = h2 * cov(z,w)
% % z = trait values before selection
% % w = relative fitness, or survival_probability / mean(surv_prob)
% if doing_DHWadaptation_covZw == 1
% %      for s = 1:nb_coral_types
% %          find_colonies =  find(coral_cm2>0);
% %          F_list = ones(length(metapop(n).coral(s).cover_cm2(find_colonies)),1); % everyone is fit
% %          [RESULT.coral_total_fecundity(n,t+1,s), fec_Pcol, inds_gravid] = f_estimate_fecundity(metapop(n).coral(s).cover_cm2(find_colonies), F_list,...
% %              CORAL.fecund_min_size(s), CORAL.fecund_a(s), CORAL.fecund_b(s));
% %      end
% %     [RESULT.coral_total_fecundity(n,t+1,s), fec_Pcol, inds_gravid] = f_estimate_fecundity(metapop(n).coral(s).cover_cm2(find_colonies), F_list,...
% %     CORAL.fecund_min_size(s), CORAL.fecund_a(s), CORAL.fecund_b(s));
% % function [fecundity, fecundity_per_colony, I] = f_estimate_fecundity (coral_cm2, F_list, fecund_min_size, a, b)
% % % Select colony sizes with 100% gravid
% % I=find(coral_cm2>= fecund_min_size);
% % adult_coral_sizes = coral_cm2(I);
% % adult_relfitness = F_list(I);
% % % Allometric relationship based on Hall and Hughes 1996:
% % % Egg volume = exp(a+b*log(size))
% % all_egg_volumes = exp(a + b*log(adult_relfitness.*adult_coral_sizes)) ; % mm3 of eggs produced by each colony
% % fecundity = floor(sum(sum(all_egg_volumes))/0.1) ; %0.1 mm3 is the average volume of an egg
% % fecundity_per_colony = floor(all_egg_volumes/0.1) ; %0.1 mm3 is the average volume of an egg
% %     HT_adult_g = metapop(n).coral(s).heat_tolerance(find_colonies(inds_gravid)); % first subset only the gravid colonies HTs,
% %     rep_fec = ceil(50*fec_Pcol/max(fec_Pcol)); % representative fecundity of each colony
% %     HT_larv = repelem(HT_adult_g, rep_fec); % then repeat the HTs of each colony by the fecundity
% %     REEF(n).coral(s).HT_pool_OUT.N(t+1)=RESULT.coral_total_fecundity(n,t+1,s); % Number of larvae produced
% %     REEF(n).coral(s).HT_pool_OUT.mu(t+1)=full(mean(HT_larv)); % average HT of larval pool
% %     [RESULT.coral_total_fecundity(n,t+1,s), fec_Pcol, inds_gravid] = f_estimate_fecundity(metapop(n).coral(s).cover_cm2(find_colonies), F_list,...
% %     CORAL.fecund_min_size(s), CORAL.fecund_a(s), CORAL.fecund_b(s));
% %     HT_adult_g = metapop(n).coral(s).heat_tolerance(find_colonies(inds_gravid)); % first subset only the gravid colonies HTs,
% %     rep_fec = ceil(50*fec_Pcol/max(fec_Pcol)); % representative fecundity of each colony
% %     HT_larv = repelem(HT_adult_g, rep_fec); % then repeat the HTs of each colony by the fecundity
% %     REEF(n).coral(s).HT_pool_OUT.N(t+1)=RESULT.coral_total_fecundity(n,t+1,s); % Number of larvae produced
% %     REEF(n).coral(s).HT_pool_OUT.mu(t+1)=full(mean(HT_larv)); % average HT of larval pool
% 
% for s = 1:nb_coral_types
%     
% end
%     z = full(heat_tolerance(coral_cm2 > fecund_min_size)); % same as HT_adult_g
%     w = 1 - full(prob_whole_mortality(coral_cm2 > fecund_min_size));
%     % weight the lists of trait values and fitness values by the fecundity
%     % of each coral colony. To do this capture corresponding a, b, and colony sizes
%     adult_coral_sizes = full(coral_cm2(coral_cm2 > fecund_min_size));
%     fecund_a = full(fecund_a(coral_cm2 > fecund_min_size));
%     fecund_b = full(fecund_b(coral_cm2 > fecund_min_size));
%     all_egg_volumes = exp(fecund_a + fecund_ab*log(1.*adult_coral_sizes)) ; % mm3 of eggs produced by each colony
%     fecundity = floor(sum(sum(all_egg_volumes))/0.1) ; %0.1 mm3 is the average volume of an egg
%     fecundity_per_colony = floor(all_egg_volumes/0.1) ; %0.1 mm3 is the average volume of an egg % same as fec_Pcol
% 
%     rep_fec = ceil(50*fecundity_per_colony/max(fecundity_per_colony)); % representative fecundity of each colony
%     HT_larv = repelem(z, rep_fec); % then repeat the HTs of each colony by the fecundity
% %     REEF(n).coral(s).HT_pool_OUT.N(t+1)=RESULT.coral_total_fecundity(n,t+1,s); % Number of larvae produced
% %     REEF(n).coral(s).HT_pool_OUT.mu(t+1)=full(mean(HT_larv)); % average HT of larval pool
% 
%     rw = w / mean(w);
% %     scatter(z,rw);
% %     cv = cov(z,rw);
% %     cv = cv(2,1);
% %     cv = sum((z-mean(z)).*(rw-mean(rw)))./length(z);
% %     R = CORAL.DHWbleaching_mortality_h2 * cv
%     R = CORAL.DHWbleaching_mortality_h2 * sum((z-mean(z)).*(rw-mean(rw)))./length(z);
% end

id_dead = id1 ;
id_dead(rand_mort1 > prob_whole_mortality) = 0 ; % exclude the survivors
coral_loss = coral_cm2.*id_dead ;

algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_loss,2) ;
% coral_cm2(id_dead==1) = - coral_cm2(id_dead==1);  % now dead (negatives)
coral_cm2(id_dead==1) = 0; 
% ht = heat_tolerance(heat_tolerance~=0); mean(ht,'all')
heat_tolerance(id_dead==1) = 0; 
% ht = heat_tolerance(heat_tolerance~=0); mean(ht,'all')
if doing_coral_age == 1
    coral_age(id_dead==1) = 0; % update so dead corals have age zero
end
total_mortality = sum(sum(id_dead))/sum(sum(id1));

%________________________________
%
% Partial mortality
%________________________________

% These are used to update the status of each coral. If a coral has never been bleached
% and then survives a bleaching event it's status changes from 0 to 1.
% NOTE from JH: the record of previous bleaching doesn't  have any effect on likelihood
% of partial mortality whereas it does on total mortality, doesn't make sense.
% NOTE from YM: we now (01/2015) apply the same reduction to the probability of partial mortality 

id2 = id1 - id_dead ;
rand_mort2 = rand(size(id2)) ; % Generates random probability for adol + adults

if doing_DHWbleaching == 1
    
    % following previous implementation
%     prob_partial_mortality = id2.*(sensitivity_bleaching .* bleaching_mortality_prob);
    % for either NO taxa sensitivities, or sensitivity defined by shift
    % heat tolerance
    prob_partial_mortality = id2.*(bleaching_mortality_prob);
    
else
    prob_partial_mortality = id2.*(1-(1-sensitivity_bleaching.*MORT).^bleaching_partial_offset);
    % prob_partial_mortality = id2.* sensitivity_bleaching .* MORT * bleaching_partial_offset;
end


id_part = id2 ;
id_part(rand_mort2 > prob_partial_mortality) = 0 ; 
bleach_extent = floor(extent_bleaching.* coral_cm2 .* id_part) ; 
algal_cm2(:,1) = algal_cm2(:,1) + sum(bleach_extent, 2) ;
coral_cm2 = coral_cm2 - bleach_extent ;

%________________________________
%
% Clade switching
%________________________________
if doing_clades == 1
    rand_switch = rand(size(id2)) ; % Generates random probability of switching to the thermally-tolerant clade (clade 2)
    clade(rand_switch < proba_switching & coral_cm2>0) = 2 ; % NOTE THIS INDEPENDENT OF BLEACHING MORTALITY
end


%%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild(coral_cm2, surface_cm2, volume_cm3, coral_age, colony_ID, clade, heat_tolerance, species_ID, nb_coral_types, doing_clades, doing_DHWbleaching, doing_coral_age, doing_3D);
% ht = coral(3).heat_tolerance(coral(3).heat_tolerance~=0); mean(ht,'all')

for a=1:size(algal_cm2, 2) 
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end

% total_coral_loss = sum(sum(coral_loss,2)) + sum(sum(bleach_extent, 2)) ;

% NEW Jul 2020: keep track of losses per species
count = 1;
total_coral_loss = zeros(1,nb_coral_types);

for s = 1:nb_coral_types
    select_sp = count:(count+species_ID(s)-1);
    total_coral_loss(s) = sum(sum(coral_loss(:,select_sp)))+sum(sum(bleach_extent(:,select_sp)));
    count = count+species_ID(s);
end    
    
