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

function [zo_mu, zo_sd] = f_additive_genetic_trait_inheritance_fast(z, f, h2)
if sum(f) >= 10000
        z_cond = full(repelem(z, ceil(f*1e3/sum(f)))); % subsample 10,000 from the egg pool - changed from round to ceil to avoid pulling zero from a colony
else
    z_cond = full(repelem(z, f));
end
x = linspace(-10,10,1e3+1);
z_std = (z_cond - mean(z_cond))/std(z_cond);
% g_pdf = pdf('Normal',x,mean(h2 * z_std),std(h2 * z_std));
% g_pdf = g_pdf/sum(g_pdf);
% g_pdf = conv(g_pdf,g_pdf,'same');
% g_pdf = g_pdf/sum(g_pdf);
% cs = cumsum(g_pdf);
% ind = find(cs >= 0.1587,1); % 1 SD below mean is cumulative probability density integral of 0.1587 = (1-0.6827)/2 for distribution which integrates to 1
% g_pm_sd = -x(ind-1) -(x(ind)-x(ind-1)).*(0.1587 - cs(ind-1)) ./ (cs(ind) - cs(ind-1));  % add the proportional distance from x(ind-1) to the true SD
% zo_mu = (x(g_pdf == max(g_pdf))/h2)*std(z_cond) + mean(z_cond); %zo_mu
% zo_sd = ((x(g_pdf == max(g_pdf)) + g_pm_sd)/h2)*std(z_cond); %zo_sd
if h2 > 0
    g_pdf = pdf('Normal',x,mean(h2 * z_std),std(h2 * z_std));
    g_pdf = g_pdf/sum(g_pdf);
    g_pdf = conv(g_pdf,g_pdf,'same');
    g_pdf = g_pdf/sum(g_pdf);
end
if h2 < 1 
    e_pdf = pdf('Normal',x,mean((1-h2)*z_std),1-h2);
    e_pdf = e_pdf/sum(e_pdf);
end
if h2 < 1 && h2 > 0
    z_std_pdf = conv(g_pdf,e_pdf,'same');
elseif h2 == 1
    z_std_pdf = g_pdf;
elseif h2 == 0
    z_std_pdf = e_pdf;
end
cs = cumsum(z_std_pdf);
ind = find(cs >= 0.1587,1); % 1 SD below mean is cumulative probability density integral of 0.1587 = (1-0.6827)/2 for distribution which integrates to 1
z_std_sd = -x(ind-1) -(x(ind)-x(ind-1)).*(0.1587 - cs(ind-1)) ./ (cs(ind) - cs(ind-1));  % add the proportional distance from x(ind-1) to the true SD
zo_mu = (x(z_std_pdf == max(z_std_pdf)))*std(z_cond) + mean(z_cond); %zo_mu
zo_sd = (z_std_sd)*std(z_cond); %zo_sd