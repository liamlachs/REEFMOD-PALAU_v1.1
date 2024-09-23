%__________________________________________________________________________
%
% REEFMOD-PALAU MAIN SCRIPT
%
% This is updated version of REFMOD-PALAU developed from REEFMOD-GBR v6.8 (Y-M Bozec) to be used 
% to capture adaptation through natural selection for heat tolerance. Implements individual
% tracking of heat tolerance and trait inheritance, and includes demographic sensitivity analyses.
%
% Updated by Liam Lachs, liamlachs@gmail.com, 02/2024
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 08/2022
%__________________________________________________________________________

function MAIN_REEFMOD_PAL(h2, ssp, gcm, ss, BH_a, BH_b)
% For testing
% h2=0; % heritability
% ssp=2; % greenhouse gas emissions scenario (Shared Socioeconoic Pathway)
% gcm=16; % Global Circulation Model
% ss = 0; % Self-seeding adjustment factor, default value is 0
% BH_a = 0; % Beverton-Holt function density dependent settlement adjustment factor, default value is 0
% BH_b = 0; % Beverton-Holt function density dependent settlement adjustment factor, default value is 0



clearvars -except h2 ssp gcm ss BH_a BH_b

SaveDir ='outputs2\';

NB_SIMULATIONS = 1; % Number of repeated runs

% NB_TIME_STEPS has to be an even number <= 26+158 (length of projected DHW time series is 79 years)
% We always run the hindcast (26 time steps) before future projections
% Initialisation = winter 2007

NB_TIME_STEPS = 26+158; % HINDCAST+FORECAST summer 2008 - winter 2099

OutputName = 'A2_FORECAST_PAL_h2'; 

options = [0 1 0 0 0]; % bleaching only - see options below

% select the Global Circulation Model for climate change projection 
GCM = gcm; % 1=ACCESS-CM2, 2=ACCESS-ESM1-5, 3=BCC-CSM2-MR, 4=CanESM5, 5=CESM2, 
         % 6=CNRM-ESM2-1, 7=EC-Earth3, 8=EC-Earth3-Veg, 9=HadGEM3-GC31-LL,
         % 10=IPSL-CM6A-LR, 11=MIROC6, 12=MPI-ESM1-2-LR, 13=MRI-ESM2-0, 
         % 14=NESM3, 15=NorESM2-LM, 16=UKESM1-0-LL
% Select the carbon emission pathway - note SSP3 for HadGEM3-GC31-LL or NESM3
SSP = ssp; % 1=SSP1-2.6, 2=SSP2-4.5, 3=SSP3-7.0, 4=SSP5-8.5

%% --------------------------------------------------------------------------------
GCM_list = ["ACCESS-CM2"; "ACCESS-ESM1-5"; "BCC-CSM2-MR"; "CanESM5"; "CESM2";
            "CNRM-ESM2-1"; "EC-Earth3"; "EC-Earth3-Veg"; "HadGEM3-GC31-LL";
            "IPSL-CM6A-LR"; "MIROC6"; "MPI-ESM1-2-LR"; "MRI-ESM2-0";
            "NESM3"; "NorESM2-LM"; "UKESM1-0-LL"];
SSP_list = ["126"; "245"; "370"; "585"];
OPTIONS.GCM = GCM_list(GCM);
OPTIONS.SSP = SSP_list(SSP);

% Stressor options: yes(1)/no(0)
OPTIONS.doing_cyclones = options(1);
OPTIONS.doing_bleaching = options(2) ; 
OPTIONS.doing_COTS = options(3);
OPTIONS.doing_WQ = options(4);

OPTIONS.SatGCM_transition = 1; % 1= 5-year weighted transition(satellite to GCM DHWs), 0 is clean break

OPTIONS.doing_size_frequency = 1; % for tracking population size structure (incl. juveniles)
OPTIONS.doing_adaptation = 0 ;
OPTIONS.adaptation_parms = [ 1 1 1.5 ]; % sigma cold/sigma hot/esd
OPTIONS.doing_restoration = options(5) ;

OPTIONS.doing_DHWbleaching = 1; %0 ;
OPTIONS.doing_DHWadaptation = 1;  %0 ;
% For the selection differential (S) in the breeders equation, Zs is the
% selected parents and Zp is the population mean. Response to selection
% (R) = h2 * (Zs - Zp). 
% (0) Zp is from the previous year (classic breeders equation)
% (1) Zp is from the ancestral population. 
OPTIONS.doing_DHWadaptation_Zp = 1; 
if OPTIONS.doing_DHWadaptation_Zp; BreedEq=''; else; BreedEq='_classicBreedEq'; end
OPTIONS.doing_coral_age = 0; 
OPTIONS.DHWbleaching_mortality_h2 = h2; % a value between 0.0 and 1.0
OPTIONS.doing_DHWbleaching_TaxaSensivities = 0;
OPTIONS.doing_Truncated = 0; % maximum heat tolerance possible is (1) the max level seen today, or (2) free to go as high
if OPTIONS.doing_Truncated; Trunc='trunc'; else; Trunc='free'; end
OPTIONS.randomize_DHW_chronology = 0; % to randomise the sequence of DHW events (DHWr) or keep in the native sequence
if OPTIONS.randomize_DHW_chronology; DHWrandom='_DHWr'; else; DHWrandom=''; end

OPTIONS.doing_1AcroporaOnly = 1 ;

% ecological tests (ET) - 
% Play with current parameters by changing them by +-X%
% (e.g., for +-20%, from -0.2*value to +0.2*value)
% multipliers can go from 0 to 1 (change value by from 0% to +- 100%)
if  ss == 0; OPTIONS.doing_selfseed_test = 0; else; OPTIONS.doing_selfseed_test = 1; end % whether to test/alter self seeding rate
if BH_a== 0; OPTIONS.doing_BH_alpha_test = 0; else; OPTIONS.doing_BH_alpha_test = 1; end % whether to test/alter density-dependent settlement rate
if BH_b== 0; OPTIONS.doing_BH_beta_test  = 0; else; OPTIONS.doing_BH_beta_test  = 1; end % whether to test/alter density-dependent settlement rate
% self-seeding (ss)
OPTIONS.coral_min_selfseed_multiplier = ss; % set at top of script (default is 0.28)
if OPTIONS.doing_selfseed_test; ET_ss=['_ETs' regexprep(regexprep(char(string(ss)),'0',''),'-','n')]; else; ET_ss=''; end
% alpha of density-dependent settlement (a)
OPTIONS.BH_alpha_multiplier = BH_a; % set at top of script 
if OPTIONS.doing_BH_alpha_test; ET_bha=['_ETa' regexprep(regexprep(char(string(BH_a)),'0',''),'-','n')]; else; ET_bha=''; end
% beta of density-dependent settlement (b)
OPTIONS.BH_beta_multiplier = BH_b; % set at top of script 
if OPTIONS.doing_BH_beta_test; ET_bhb=['_ETb' regexprep(regexprep(char(string(BH_b)),'0',''),'-','n')]; else; ET_bhb=''; end
% edit naming
ET_addon = regexprep([ET_ss ET_bha ET_bhb], '_ET','');
if ~isempty(ET_addon); ET_addon = ['_ET' ET_addon]; end
% ET_addon

% Additional analyses added upon manuscript revision
% NOTE only have 1 of the below three tests turned on at any one time
% (1) Test the recovery rate (RR) after a catastrophic disturbance
OPTIONS.doing_RR_test = 0; % set as 0 (off) or 1 (on)
if OPTIONS.doing_RR_test
    % 8C-weeks to look at decline (vs 1998)... DHW was 5.8 in 1998 and 6.8 in 2010
    % 16 C-weeks to look at recover from catastrophic disturbance
    OPTIONS.RR_test_DHWpert = 8; 
    RR=['_RecRate' char(string(OPTIONS.RR_test_DHWpert))];
    NB_TIME_STEPS = 30+30; % 15 years burn in and 15 years of recovery
else
    OPTIONS.RR_test_DHWpert = NaN;
    RR='';
end
% (2a) what if bleaching mortality is reduced by acclimatisation from pre-pulses?
OPTIONS.doing_PPacclim = 0; % set as 0 (off) or 1 (on)
if OPTIONS.doing_PPacclim; PPacclim='_PPacclim'; else; PPacclim=''; end
% (2b) what if thermal tolerance (TT) can change through acclimatisation (acclim)?
OPTIONS.doing_TTacclim = 1; % set as 0 (off) or 1 (on)
if OPTIONS.doing_TTacclim; TTacclim='_TTacclim'; else; TTacclim=''; end
% (3) what if there is was demographic rescue (DR) from outside Palau with
% incoming of same HT as resident Palau-wide mean (assumption that outside populations have had similar evolutionary process, but better demographic status (e.g. population size)
OPTIONS.doing_DemResc = 0; % set as 0 (off) or 1 (on)
if OPTIONS.doing_DemResc; DemResc='_DemRescX2'; else; DemResc=''; end
% end of ecological tests

% Below options are used to force simulations with specific starting conditions (eg, building LUT for the RRAP-RE)
% Keep them empty if of no use
OPTIONS.init_coral_cover = []; %0.01*ones(1,6); % as proportional cover (vector of 6 values)
OPTIONS.init_sand_cover = []; % as proportional cover 
OPTIONS.init_rubble_cover = []; % as proportional cover 
OPTIONS.ssc = []; %0.1; %in mg/L

%% CoTS control
OPTIONS.doing_COTS_control= 0; % Note control is set to start in 2019 (after 23 time steps)
OPTIONS.CoTS_control_scenarios = csvread('parsList2.csv', 0, 1);
% Caro: 5 boats, whole GBR, start control in 2019. New list adjusted for matching Coconet
% Option names in parsList2.csv. Second value is spatial strategy (14: whole GBR under CoTS control).
% Fourth value is number of boats (5)
OPTIONS.CoTS_control_scenarios(4)=5;

%% Outplanting
% Set the restoration effort: maximimum number of reefs where coral outplanting is undertaken at random at each time step
RESTORATION.nb_reefs_outplanted = Inf ;
% Set to Inf if unlimited OR if specific reefs are restored (listed in SETTINGS_RESTORATION
% If 0, outplanting cannot happen. For the counterfactual, set to Inf with outplanted_density = 0 for ghost deployment

RESTORATION.total_nb_outplants = Inf; % Max number of outplants available for all reefs at each time step.
% Set to Inf if outplant density is the driver

RESTORATION.outplanted_density = 0;  %only in the case of fixed density of outplants (ignores RESTORATION.total_nb_outplants)

RESTORATION.doing_coral_outplanting = zeros(1,NB_TIME_STEPS); % if zero, no coral deployment at time step for all reefs
RESTORATION.doing_coral_outplanting(1,38:2:46) = 1 ; % set to 1 to indicate outplanting: first in summer 2026 and last in summer 2030 (BCA)
% RESTORATION.doing_coral_outplanting(1,[38:2:46 78:2:86 118:2:126 158:2:166] ) = 1 ; % outplanting starts 2026-2030, 2046-2050, 2066-2070, 2086-2090

RESTORATION.thermal_tolerance_outplants = 0 ; % DegC above thermal tolerance of native corals (WITH GENETICS)
RESTORATION.DHW_tolerance_outplants = 4 ; % DHW tolerance relative to native corals (WITHOUT GENETICS)
RESTORATION.tradeoff_growth_outplant = 1; %0.6; % 40% reduction of growth rate relative to native corals

RESTORATION.outplant_species = [1 2 3 4 5 6]; % ID of outplanted coral types - this will create as many 'new' types after the 6 default types
RESTORATION.outplant_species_prop = [0.02 0.14 0.14 0 0.7 0]; % Taxonomic composition of outplanted corals: must sum to 1
RESTORATION.outplant_diameter_mean = [2.56 2.56 2.56 1.41 1.41 1.41] ; % mean diameter (in cm) of outplants of outplants of each deployed type
RESTORATION.outplant_diameter_sd = [0.26 0.26 0.26 0.14 0.14 0.14] ; % sd diameter (in cm) of outplants of outplants of each deployed type

%% Rubble stabilisation
% Set the restoration effort: number of reefs where rubble is stabilised at each time step
RESTORATION.nb_reefs_stabilised = 0 ; % (if 0 rubble stabilisation cannot happen) 
% Set the timing of intervention (if 1 intervention is deployed at step t, if 0 no intervention at t)
RESTORATION.doing_rubble_stabilisation = zeros(1,NB_TIME_STEPS);
% RESTORATION.doing_rubble_stabilisation(1,38:1:end) = 1 ; % set to 1 to do outplanting: first in summer 2026 and last in summer 2030

%% Larval enrichment
RESTORATION.nb_reefs_enriched = Inf ; % max number of reefs where larval enrichment is undertaken at each time step
% Set to Inf if unlimited OR if specific reefs are restored (listed in SETTINGS_RESTORATION)
% If 0, enrichment cannot happen. For the counterfactual, set to Inf with total_nb_larvae = 0 for ghost deployment

RESTORATION.total_nb_larvae = 1e6; % Max number of 'larvae' (ie, 1 yr old corals) available at each time step. Set to Inf if unlimited.

RESTORATION.doing_larval_enrichment = zeros(1,NB_TIME_STEPS);
% RESTORATION.doing_larval_enrichment(1,38:2:46) = 1 ; % set to 1 to do outplanting

RESTORATION.DHW_tolerance_larvae = 0;

%% SRM (cooling)
% Cloud brightening from RRAP feasibility study (needs revision)
RESTORATION.doing_cooling = 0 ;
RESTORATION.cooling_factor = 0 ; % otherwise [-0.3 ; -0.7 ; -1.3];

% Fogging for 2022 intervention simulations
RESTORATION.nb_reefs_fogged = Inf;
% Select reefs for the deployment of fogging: 695-Moore; 697-Elford; 698-Briggs ; 969-Milln; 970-Thetford
RESTORATION.fogged_reef_ID = [969 970]; % 1 fogging unit (4.5 km2) covering 4.6 km2 (equivalent 2D reef areas)
% RESTORATION.fogged_reef_ID = 695; % 2 fogging units (9 km2) covering 8.7 km2 (equivalent 2D reef areas)
% RESTORATION.fogged_reef_ID = [695 697 969 970]; % 5 fogging units (22.5 km2) covering 23 km2 (equivalent 2D reef areas)

RESTORATION.doing_fogging = zeros(1,NB_TIME_STEPS); %if zero don't do fogging at time step
% RESTORATION.doing_fogging(1,1:2:end) = 1 ; % set to 1 to do fogging - only in summer!!
RESTORATION.bleaching_mortality_under_fogging = 0.8 ; % fogging reduces bleaching mortality by 20%

%% Saving options
if NB_TIME_STEPS > 30   
    if OPTIONS.doing_DHWadaptation == 0
        OPTIONS.OutputFileName = [SaveDir OutputName '_' char(OPTIONS.GCM) '_' char(OPTIONS.SSP) '.mat'];
    else
        if OPTIONS.doing_RR_test
            OPTIONS.OutputFileName = [SaveDir OutputName sprintf('%1.1f',h2) BreedEq ET_addon RR '_' Trunc '.mat'];            
        else
            OPTIONS.OutputFileName = [SaveDir OutputName sprintf('%1.1f',h2) BreedEq ET_addon DemResc TTacclim PPacclim '_' Trunc '_' char(OPTIONS.GCM) '_' char(OPTIONS.SSP) DHWrandom '.mat'];
        end
    end
else
    OPTIONS.OutputFileName = [SaveDir OutputName '.mat'];
end

%% --------------------------------------------------------------------------------
OUTPUTS = struct('REEF', [],'RESULT', [],'RECORD', []);
TEMP_META = struct('META', []);
run_id = 1;



% parfor run_id = 1:NB_SIMULATIONS

for run_id = 1:NB_SIMULATIONS
%     run_id = 1
    run_id

    [meta, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, NB_TIME_STEPS, run_id);
      
    OUTPUTS(run_id).RESULT = RESULT ;
    OUTPUTS(run_id).RECORD = RECORD ;
    OUTPUTS(run_id).REEF = REEF ;
    TEMP_META(run_id).META = meta ;
    
end

META = TEMP_META(1).META; % Keep only one META because common to all simulations
clear TEMP_META ADAPT run_id meta RECORD REEF RESULT GCM GCM_list options SaveDir OutputName SSP SSP_list


%% --------------------------------------------------------------------------------
%% Memory allocation
% 1) coral outputs
coral_cover_per_taxa = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
coral_larval_supply = coral_cover_per_taxa;
nb_coral_offspring = coral_cover_per_taxa;
nb_coral_recruit = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'uint16');
    
if OPTIONS.doing_size_frequency == 1   
    nb_coral_juv = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 2, 'uint16') ;
    nb_coral_adol = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 3, 'uint16') ;
    nb_coral_adult = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 6, 'uint16') ; 
end

if OPTIONS.doing_DHWbleaching == 1   
    nb_super_tolerant_colonies = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single');
    nb_fecund_colonies = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single');
    heat_tolerances_mu_fecund_colonies = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single');
    heat_tolerance_mu = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single') ;
    heat_tolerance_sd = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single') ;
    heat_tolerance_975 = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 'single') ;
end

% 2) Coral cover loss following stressors
coral_cover_lost_bleaching = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, META.nb_coral_types, 'single');
coral_cover_lost_cyclones = coral_cover_lost_bleaching;
coral_cover_lost_COTS = coral_cover_lost_bleaching;

% 3) Stress records
record_applied_DHWs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
record_applied_bleaching_mortality  = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
record_applied_cyclones = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'uint16');

% 4) Restoration records
if OPTIONS.doing_restoration==1
    
    % Total nb of outplants per reef per time step
    record_total_outplants_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, length(META.outplant_species),'uint16');
    record_outplanted_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    record_total_larvae_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, length(META.enriched_species),'uint16');
    record_enriched_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');  
    
    record_rubble_pct2D_stabilised = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
    record_stabilised_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    coral_cover_per_taxa_restored_sites = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
end

% 5) Other variables
rubble = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
nongrazable = zeros(NB_SIMULATIONS, META.nb_reefs,'single');
macroTurf = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroEncrustFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroUprightFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');

% 6) CoTS outputs
if OPTIONS.doing_COTS == 1
    COTS_mantatow = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'single');
    COTS_densities0 = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 16, 'single'); % 16 age classes
    COTS_settler_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'uint16');
    COTS_larval_supply = COTS_mantatow;
    COTS_larval_output = COTS_mantatow;
end

if OPTIONS.doing_COTS_control == 1 && META.nb_time_steps > 23
    COTS_CONTROL_culled_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % gives which reefs were culled or not (1/0)
    COTS_CONTROL_remaining_dives = zeros(NB_SIMULATIONS, META.nb_time_steps); %remaining number of control dives available after intervention
    COTS_CONTROL_culled_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % density extracted by the intervention (adults)
end

%% Populate outputs
for simul = 1:NB_SIMULATIONS
    
    coral_cover_per_taxa(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D));
    coral_larval_supply(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_larval_supply)); % nb of incoming larvae per unit of reef area (400m2)
    nb_coral_recruit(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_settler_count));
    nb_coral_offspring(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_total_fecundity)); % nb of larvae produced per unit of reef area (400m2)
    
    if OPTIONS.doing_size_frequency == 1
        nb_coral_juv(simul,:,:,:,:)= squeeze(cat(4,OUTPUTS(simul).RESULT.coral_juv_count(:,:,:,:)));
        nb_coral_adol(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adol_count(:,:,:,:)));
        nb_coral_adult(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adult_count(:,:,:,:)));
    end
    
    if OPTIONS.doing_DHWbleaching == 1   
        nb_super_tolerant_colonies(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.nb_super_tolerant_colonies)) ;
        nb_fecund_colonies(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.nb_fecund_colonies)) ;
        heat_tolerances_mu_fecund_colonies(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.heat_tolerances_mu_fecund_colonies)) ;
        heat_tolerance_mu(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.heat_tolerances_mu)) ;
        heat_tolerance_sd(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.heat_tolerances_sd));
        heat_tolerance_975(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.heat_tolerances_975));
    end
    
    coral_cover_lost_bleaching(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_bleaching);
    coral_cover_lost_cyclones(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_cyclones);
    coral_cover_lost_COTS(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_COTS);
    
    record_applied_DHWs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_DHWs);
    record_applied_bleaching_mortality(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_bleaching_mortality);
    record_applied_cyclones(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.hurricane_events);
    
    if  OPTIONS.doing_restoration==1
        record_total_outplants_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_outplanted);
        record_outplanted_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.outplanted_reefs);
        record_total_larvae_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_enriched);
        record_enriched_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.enriched_reefs);
        record_rubble_pct2D_stabilised(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.rubble_cover_pct2D_stabilised(:,1:end));
        record_stabilised_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.stabilised_reefs);
        
        coral_cover_per_taxa_restored_sites(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D_restored_sites));
    end
    
    rubble(simul,:,:)= squeeze(OUTPUTS(simul).RESULT.rubble_cover_pct2D);
    nongrazable(simul,:) = squeeze(cat(4,OUTPUTS(simul).REEF.nongrazable_substratum));
    
    macroTurf(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,4)));
    macroEncrustFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,2)));
    macroUprightFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,3)));
      
    if OPTIONS.doing_COTS == 1
        % Assuming 0.6 CoTS per grid ~ 0.22 CoTS per tow
        % (0.22 per tow is equivalent to 1500 COTS per km2 (Moran & De'ath 92), so that 1 COTS per grid (400m2) is equivalent to 0.22*2500/1500
        COTS_mantatow(simul,:,:) = (0.22/0.6)*squeeze(cat(4,OUTPUTS(simul).RESULT.COTS_total_perceived_density));
        COTS_densities0(simul,:,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities); % Density for 400m2
        COTS_settler_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_settler_densities); % Density for 400m2
        COTS_larval_supply(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_supply); % Density for 400m2
        COTS_larval_output(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_output); % Density for 400m2
    end
    
    if OPTIONS.doing_COTS_control == 1 && META.nb_time_steps > 23
        COTS_CONTROL_culled_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_reefs);
        COTS_CONTROL_remaining_dives(simul,:) = squeeze(OUTPUTS(simul).RESULT.COTS_control_remaining_dives);
        COTS_CONTROL_culled_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_density);
    end

end 

% New (08/2021): only record COTS densities by yearly classes to reduce output size
% COTS_densities = COTS_densities0(:,:,:,1:2:end)+COTS_densities0(:,:,:,2:2:end);
% 09/2021: now just sum across all juveniles and across all adults
if OPTIONS.doing_COTS == 1
    COTS_juv_densities = sum(COTS_densities0(:,:,:,1:(META.COTS_adult_min_age-1)),4);
    COTS_adult_densities = sum(COTS_densities0(:,:,:,META.COTS_adult_min_age:end),4);
end

clear OUTPUTS simul s COTS_densities0 NB_TIME_STEPS

save (OPTIONS.OutputFileName)
% end