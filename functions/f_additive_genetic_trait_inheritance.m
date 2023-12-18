% REEFMOD-PAL script to compute the offspring trait distribution on a reef
% given the adult parental distribution
%
% Liam Lachs, l.lachs2@newcastle.ac.uk, 09/2022
%
%     z: trait values of the adults in the population after selection
%     f: fecundity of each individual adult
%    h2: narrow-sense heritability of the trait
%        (CORAL.DHWbleaching_mortality_h2)
% zo_mu: mean of the offspring trait distribution
% zo_sd: SD of the offspring trait distribution
%__________________________________________________________________________

function [zo_mu, zo_sd] = f_additive_genetic_trait_inheritance(z, f, h2)


%% Doing a proper job of inhiting heat tolerances into the larval pool
%% 0) Create a conditional trait distribution
% For corals this implies making an egg pool based on colony fecundity,
% where each egg inherits exactly it's mothers trait
% The adult trait distribution weighted by fecundity
% a random sample of this is fine otherwise the vecotr is too big
z_cond = full(repelem(z, f));
if length(z_cond) >= 10000
    z_cond = datasample(z_cond, 10000); % subsample 10,000 from the egg pool
end

%% 1) Standadise the parental phenotype (z) post selection
% this will allow us to then use a PO-style regression to seperate out just
% the environmental component
% HT_adult_ps = normrnd(6.86, 1, 10000,1);
z_mu = mean(z_cond);
z_sd = std(z_cond);
z_std = (z_cond - z_mu)/z_sd;
clear z_cond
% mean(z_std);

%% 2) separate out the genetic component of the phenotypic distribution
% z_std_offspring = B0 + B1 * z_std_parent;
% where z_std_offspring can be replaced by the genotypic component
% B0 = 0, and B1 = narrow-sense heritability (h2)
% as such g = h2 * z_std
g = h2 * z_std;
% clear z_std

%% 3) parent midpoint (pm) genotype distibution
% this is achieved by doing a convolution of the genotype probability
% density function on itself. Following Coulson. Maybe this can also be
% done with a non-normal pdf?
x = linspace(-10,10,1e3+1);
if h2 > 0
    g_pdf = pdf('Normal',x,mean(g),std(g));
    g_pdf = g_pdf/sum(g_pdf);
    g_pdf_pm = conv(g_pdf,g_pdf,'same');
    % clear g_pdf
    % sum(g_pdf_pm) % should integrate to 1
    % figure
    % scatter(x,g_pdf); hold on
    % scatter(x,g_pdf_pm)
    % x(find(g_pdf_pm == max(g_pdf_pm))) mean should be zero
    
    % mean
    g_pm_mu = x(g_pdf_pm == max(g_pdf_pm)); % round(mean(x(randsample(length(x),1e6,true,g_pdf_pm))),2);
    % SD - 1 million random samples give SD precision to 2 decimal places
    g_pm_sd = round(std(x(randsample(length(x),1e6,true,g_pdf_pm))),2);
end
clear x

%% 4) add segregation variance
% variance due to the random distribution of either copy (dominant Y or
% recessive y) of a homozygous gene from parent to gamete
% following Coulson this is achieved by convolving th parent midpoint
% genotype distribution with a normal PDF of mean 0 and sd ?
% not implementing as unsure of effect size

% %% 5a) Add back the environmental variance component to get the offspring phenotype (standardised)
% % this is achieved in the opposite way to separating out the additive
% % genetic component (g=h2*z_std), instead we use z_std = g/h2
% % zo_std = normrnd(g_pm_mu,g_pm_sd,1e4,1)/h2;
%
%
% %% 6a) back-standardise to get the true offspring phenotypic trait distribution
% % back transform based on the original standardisation of parent trait
% % values (z_std = (z - z_mu)/z_sd)
% % z = z_std * z_sd + z_mu
% % zo = zo_std * z_sd + z_mu;
% % mean(zo)
% % std(zo)
%
% % or I can skip the random number generation in step 5 and instead use the
% % following equations to calculate the new mean and SD
% % zo_mu = (g_pm_mu/h2)*z_sd + z_mu; %zo_mu
% % zo_sd = ((g_pm_mu + g_pm_sd)/h2)*z_sd; %zo_sd


%% 5b) Add back the environmnta varianc componnt to get the offspring phenotype
% this could also be done by convolving the g_pm distribution with the
% environmental component
x = linspace(-10,10,1e3+1);
if h2 < 1
    e = (1-h2) * z_std;
    e_pdf = pdf('Normal',x,mean(e),1-h2);
    e_pdf = e_pdf/sum(e_pdf);
end
if h2 == 1
    z_std_pdf = g_pdf_pm;
elseif h2 == 0
    z_std_pdf = e_pdf;
else
    z_std_pdf = conv(g_pdf_pm,e_pdf,'same');
end
clear e_pdf

%% 6b) back-standardise to get the true offspring phenotypic trait distribution
% back transform based on the original standardisation of parent trait
% values (z_std = (z - z_mu)/z_sd)
% z = z_std * z_sd + z_mu
% zo = zo_std * z_sd + z_mu;
% mean
z_std_mu = x(z_std_pdf == max(z_std_pdf)); % round(mean(x(randsample(length(x),1e6,true,g_pdf_pm))),2);
% SD - 1 million random samples give SD precision to 2 decimal places
z_std_sd = round(std(x(randsample(length(x),1e6,true,z_std_pdf))),2);


zo_mu = (z_std_mu)*z_sd + z_mu; %zo_mu
zo_sd = (z_std_sd)*z_sd; %zo_sd

