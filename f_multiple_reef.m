% REEFMOD-PALAU model run script
%
% Updated by Liam Lachs, liamlachs@gmail.com, 12/2023
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 08/2022
%
%__________________________________________________________________________

function [META, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, nb_time_steps, simul)

if OPTIONS.doing_1AcroporaOnly == 0
    PARAMETERS
else
    PARAMETERS_1AcroporaOnly
end
META.nb_time_steps = nb_time_steps;

%% Reef areas
load('PAL_REEF_POLYGONS_2022.mat');
PAL_REEFS.Reference_Area_km2 = ones(95,1);
 
%% Reef selection
META.reef_ID = (1:95)'; % Entire Palau

META.nb_reefs = length(META.reef_ID);
META.outside_reef_ID = [];  %note this is empty when all reefs are included
META.reef_lat = PAL_REEFS.LAT(META.reef_ID); % needed for CoTS control
META.reef_lon = PAL_REEFS.LON(META.reef_ID); % needed for CoTS control

% Define which habitat area to use
META.area_habitat = PAL_REEFS.Reference_Area_km2(META.reef_ID);

% Enlarge reef grids to deploy in specific sites (set to all reefs)
% Future version will allow setting a specific grid size to each reef individually
% META.grid_x_count = 43; % number of grid cells along the x-edge
% META.grid_y_count = 43; % number of grid cells along the y-edge
% For the BCA:
% META.grid_x_count = 35; % number of grid cells along the x-edge
% META.grid_y_count = 35; % number of grid cells along the y-edge

%% Disturbances
% NOTE: don't refine COTS parameters here because they will be erased by settings_COTS
META.doing_bleaching = OPTIONS.doing_bleaching ;
META.deterministic_bleaching = 1; % option for generating deterministic (1) or random (0) mortalities from DHWs
META.DHW_threshold = 3 ; % DHW threshold that triggers bleaching % ignored for DHWadaptation

% 5-yr weighted transition from satellite to GCM
META.SatGCM_transition = OPTIONS.SatGCM_transition;

% Adjust bleaching to align with the shallow (-2m) mortality recorded by Hughes et al. (2018) or deep (-7m) as in Baird et al. (2018)
CORAL.bleaching_depth = 1; % shallow bleaching as simulated for hindcast (Bozec et al. 2022)
% CORAL.bleaching_depth = 0.5; % deep bleaching (Baird et al. 2018 MEPS)

% set heritability as per options
CORAL.DHWbleaching_mortality_h2 = OPTIONS.DHWbleaching_mortality_h2;

META.doing_hurricanes = OPTIONS.doing_cyclones ;
META.random_hurricanes = 0; % put 0 to apply a specific scenario (1 imposes random occurence)
META.deterministic_hurricane_mortality = 0; % option for generating deterministic (1) or random (0) mortalities from cyclone cat

META.doing_COTS = OPTIONS.doing_COTS ;
META.doing_COTS_control = OPTIONS.doing_COTS_control ;
META.report_COTS_control = 0 ; % Put 1 to record the detailed results of CoTS control (will create RESULT.COTS_records)
META.randomize_initial_COTS_densities = 2 ; % Distri for random generation of COTS densities (based on CoCoNet averages)
% 1 for Gaussian, 2 for Poisson

% Define the COTS initiation box
META.reef_COTS_INIT_BOX = []; % not relevant to Palau?
% META.reef_COTS_INIT_BOX = PAL_REEFS.Reef_ID(GBR_REEFS.LAT >= -17 & GBR_REEFS.LAT <= -14.6); 
% Pratchett et al. 2014: bounds of the initiation box 

REEF_COTS =[]; % (required, even if doing CoTS)

% WQ options
META.doing_water_quality = OPTIONS.doing_WQ;
META.randomize_WQ_chronology = 0;
META.doing_Chl_forcing = 1; % with (1) or without (0) CoTS larval survival driven by Chlorophyll concentration

% Track colony size distributions
META.doing_size_frequency = OPTIONS.doing_size_frequency; % Track coral colonies sizes (SLOW!)

%% Connectivity
META.doing_coral_connectivity = 1 ;
% META.coral_immigration = 1e6*ones(1,META.nb_time_steps) ; % Forced larval input for a 400m2 area - works only if connectivity is OFF
% Simulating the 4 reefs of the Moore Reef cluster only:
% ES=ES_6mo';
% IMMI = zeros(META.nb_reefs,size(ES,2),size(ES,1));
% 
% for n=1:META.nb_reefs
%     IMMI(n,:,:) = ES';
% end
% META.coral_immigration = IMMI;

META.recruitment_type = 1; % turn into 0 for fixed recruitment (but then connect and genetics won't work)

META.doing_selfseed_test = OPTIONS.doing_selfseed_test;
META.doing_BH_alpha_test = OPTIONS.doing_BH_alpha_test;
META.doing_BH_beta_test = OPTIONS.doing_BH_beta_test;
META.doing_RR_test = OPTIONS.doing_RR_test;
META.RR_test_DHWpert = OPTIONS.RR_test_DHWpert;
META.doing_TTacclim = OPTIONS.doing_TTacclim;
META.doing_DemResc = OPTIONS.doing_DemResc;

% Parameter a of the B-H function (same for all reefs), calibrated with
% META.coral_min_selfseed = 0.28 (Bozec et al. 2022)
CORAL.BH_alpha = 15*CORAL.prop_settlers; % per m2
CORAL.BH_beta = 5*1e6*ones(6,1); % for a 400m2 reef

% Force self-seeding of coral larvae
META.coral_min_selfseed = 0.28 ; % relative proportion of produced larvae that stay on the reef (Helix experiment)

% Ecological test: Liam
if META.doing_selfseed_test == 1
    META.coral_min_selfseed = META.coral_min_selfseed * (1 + OPTIONS.coral_min_selfseed_multiplier);
end
if META.doing_BH_alpha_test == 1
    CORAL.BH_alpha = CORAL.BH_alpha * (1 + OPTIONS.BH_alpha_multiplier);
end
if META.doing_BH_beta_test == 1
    CORAL.BH_beta = CORAL.BH_beta * (1 + OPTIONS.BH_beta_multiplier);
end

%% Rubble (standard)
META.tracking_rubble = 1;
META.rubble_decay_rate = 0.128 ; % 2/3 stabilised after 4 years
META.convert_rubble = 1; % Conversion factor from coral loss to rubble cover
META.convert_rubble_lag = 6; % 3 years delay for converting coral loss due to bleaching and CoTS (delayed structural loss).

%% Reef with a rubble problem
% META.tracking_rubble = 1;
% META.rubble_decay_rate = 0.07 ;
% META.convert_rubble = 2.5;
% META.convert_rubble_lag = 2; % 1 year delay for converting coral loss due to bleaching and CoTS (delayed structural loss).


%% Grazing
REEF.herbivory = 1; % full grazing = 1
ALGAL.nb_step_algal_dynamics = 1 ; %%%%% ONLY TO SPEED-UP THE CODE WHEN FULL GRAZING (otherwise set to 6)

%% Restoration
META.doing_restoration = OPTIONS.doing_restoration;

%% DHW bleaching and adaptation (Liam)
META.doing_DHWbleaching = OPTIONS.doing_DHWbleaching;
META.doing_DHWadaptation = OPTIONS.doing_DHWadaptation ;
META.doing_DHWadaptation_Zp = OPTIONS.doing_DHWadaptation_Zp;
META.doing_coral_age = OPTIONS.doing_coral_age; % for overlapping generations in breeding equation
META.doing_DHWbleaching_TaxaSensivities = OPTIONS.doing_DHWbleaching_TaxaSensivities ;
META.doing_Truncated = OPTIONS.doing_Truncated;
META.doing_1AcroporaOnly = OPTIONS.doing_1AcroporaOnly;
META.randomize_DHW_chronology = OPTIONS.randomize_DHW_chronology;
META.doing_PPacclim = OPTIONS.doing_PPacclim; 

%% Genetic adaptation
META.doing_genetics = OPTIONS.doing_adaptation ;

if META.doing_genetics==1
    
    settings_GENETICS;
    
    META.genetics.SIGMA_COLD = OPTIONS.adaptation_parms(1) + [ 0 0 0 0 0 0 ]; % for the cold side (when temp<Topt)
    META.genetics.SIGMA_HOT = OPTIONS.adaptation_parms(2) + [ 0 0 0 0 0 0 ]; % for the hot side (when temp>Topt)
    META.genetics.esd = OPTIONS.adaptation_parms(3) + [ 0 0 0 0 0 0 ]; % SD of environmental effect on fitness (mean=0), on top of genetics.
    META.genetics.enhanced_tolerance = OPTIONS.thermal_tolerance_outplants ;
    
    % Load the pre-adapted pool of QTL
    load(['QTL_pool_sigma_c' num2str(OPTIONS.adaptation_parms(1)) '_h' num2str(OPTIONS.adaptation_parms(2))...
        '_esd' num2str(OPTIONS.adaptation_parms(3)) '.mat']);
    
    CORAL.growth_rate = CORAL.growth_rate/mean_fitness; % average fitness across the GBR after burn-in
    % this allows adjusting fitness so that mean individual growth rate is close to empirical values from literature
    % Note this average fitness is relative to the local environment corals are adapted to, ie we are not 
    % modelling latitudinal differences in growth rates (100% fitness in the far South gives the same growth rate than 
    % in the far North, while in reality coral growth declines at with latitude
    
end

%% INITIALISATION
% rng('shuffle')
rng(simul); % to get a repeatable scheme of random number generation in RAND, RANDI, RANDN

INITIALISATION

settings_PAL

if META.doing_restoration == 1
    
    settings_RESTORATION;
    
    % Then generate priority lists for each technique = list of reef ID sorted from highest to lowest priority
    % Note the list is set only once and remains the same throughout the simulation
    MY_REEFS = GBR_REEFS(META.reef_ID,:);
%     META.priority_list_Outplant = f_generate_priority_list_NEW(META.priority_option_Outplant, MY_REEFS, CONNECT);
%     META.priority_list_RubbleStab = f_generate_priority_list_NEW(META.priority_option_RubbleStab, META.reef_ID, MY_REEFS, CONNECT);
%     META.priority_list_LarvalEnrich  = f_generate_priority_list_NEW(META.priority_option_LarvalEnrich, MY_REEFS, CONNECT);
%     META.priority_list_Fogging = f_generate_priority_list_NEW(META.priority_option_Fogging, MY_REEFS, CONNECT);

    % ONLY FOR THE BCA: focus on the two reef clusters and deploy fixed density of outplants
    % Deployment areas are set in settings_RESTORATION
    META.priority_list_Outplant = find(META.coral_deployment.DeploymentArea_km2>0);
    META.priority_list_LarvalEnrich = find(META.coral_deployment.DeploymentArea_km2>0);
    META.priority_list_RubbleStab = 1:length(META.reef_ID); % do it everywhere for the moment
    
    reef_fogging_list = ismember(META.coral_deployment.Reef_ID,META.fogged_reef_ID);
    META.priority_list_Fogging = find(reef_fogging_list == 1);
end

% Adjustments in case of coral outplanting
if META.nb_coral_types > 6
    
    % Add initial cover for outplants (0%)
    init_coral_cover = [init_coral_cover zeros(META.nb_reefs, META.nb_coral_types-6)];
    
    if META.doing_COTS == 1
        % Extend the vector of CoTS preferences
        X = META.COTS_feeding_prefs;
        X = [X ; X];
        META.COTS_feeding_prefs = X/sum(X); % feeding prefs sum to 1
%         META.COTS_pref_corals = [META.COTS_pref_corals META.COTS_pref_corals(META.outplanted_species)+5]; % don't need to change this
    end
    
end

%% REFINE IF SPECIFIC FORCING CONDITIONS (eg, for building LUT for the RRAP-RE)
% (only works with single reef simulations)
if isempty(OPTIONS.init_coral_cover)==0
    init_coral_cover = OPTIONS.init_coral_cover;  
end

if isempty(OPTIONS.init_sand_cover)==0
    init_sand_cover = OPTIONS.init_sand_cover;  
end

if isempty(OPTIONS.init_rubble_cover)==0
    init_rubble_cover = OPTIONS.init_rubble_cover;  
end    
    
if isempty(OPTIONS.ssc)==0 
    CORAL_recruit_survival = (1 - 1.88*0.001*OPTIONS.ssc)^(180/40);
    CORAL_juvenile_growth = 1 - 0.176*log(OPTIONS.ssc+1);
    
    FERT = exp(4.579 - 0.010*OPTIONS.ssc)/exp(4.579); % fertilization success
    SETT = (99.571 - 10.637*log(OPTIONS.ssc+1))/99.571; % settlement success
    CORAL_larvae_production =  FERT*SETT;
    
    % Store in REEF_POP (same value for all years)
    for y=1:size(REEF_POP,2)
        
        REEF_POP(y).CORAL_recruit_survival = CORAL_recruit_survival;
        REEF_POP(y).CORAL_juvenile_growth = CORAL_juvenile_growth;
        REEF_POP(y).CORAL_larvae_production = CORAL_larvae_production;    
    end
end

MULTIPLE_REEF_SETUP

% Randomize the bleaching chronology (same randomised chronology across all
% reefs
    if META.doing_bleaching == 1 && META.randomize_DHW_chronology == 1
        % Add stochasticity to the DHW forcing from the single GCM run by
        % rearranging the sequence of DHW events within each decade (on a site-by-site basis)
        inds = 5:2:META.nb_time_steps; % index of each summer starting in 2010 (index t=5)
        inds(1:10)  = inds(randsample(1:10,10)); % For first decade randomly shuffle indexes
        inds(11:20) = inds(10+randsample(1:10,10));
        inds(21:30) = inds(20+randsample(1:10,10));
        inds(31:40) = inds(30+randsample(1:10,10));
        inds(41:50) = inds(40+randsample(1:10,10));
        inds(51:60) = inds(50+randsample(1:10,10));
        inds(61:70) = inds(60+randsample(1:10,10));
        inds(71:80) = inds(70+randsample(1:10,10));
        inds(81:90) = inds(80+randsample(1:10,10));
        for n = 1:META.nb_reefs
            REEF(n).predicted_DHWs(5:2:META.nb_time_steps) = REEF(n).predicted_DHWs(inds);
        end
    end
        

%% REFINE BACKGROUND COTS DENSITY ON SPECIFIC REEFS
if META.doing_COTS == 1 % LABL
    for n=1:length(META.reef_COTS_INIT_BOX)
        id_reef = find(META.reef_ID==META.reef_COTS_INIT_BOX(n));
        
        if isempty(id_reef) == 0
            REEF(id_reef).COTS_background_density = META.COTS_background_density_INIT_BOX;
        end
    end
end % LABL

%% APPLY CHLOROPHYLL FORCING OF COTS LARVAL SURVIVAL ONLY ON INSHORE REEFS (ONLY WORKS IF doing_Chl_forcing=1)
if META.doing_Chl_forcing == 1 && META.doing_COTS == 1
    for n=1:META.nb_reefs
        if GBR_REEFS.Shelf_position(META.reef_ID(n)) > 1
            for yr=1:size(REEF_POP,2)
                REEF_POP(yr).COTS_larvae_survival(n) = META.COTS_min_larval_survival;
            end
        end
    end
end

%% If doing COTS control 
if META.doing_COTS_control == 1
    [ META ] = f_settings_COTS_control(OPTIONS.CoTS_control_scenarios, META);
end
             
clearvars -except META REEF CORAL ALGAL CONNECT REEF_POP REEF_COTS

[RESULT, RECORD] = f_runmodel(META, REEF, CORAL, ALGAL, CONNECT, REEF_POP, REEF_COTS) ;
