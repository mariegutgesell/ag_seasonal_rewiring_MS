Data and code to reproduce analysis in Gutgesell, M.K., Guzzo, M.M, and McCann, K.S. (2025). Agriculture land-use change seasonally rewires stream food webs: a case study from headwater streams in the Lake Erie watershed. FACETS. 

Code Files:
seasonal_stream_morphology.R <-- code to calculate seasonal mean and sd of water depth, wetted width and hydraulic head
summer_nutrient_algae_mean.R <-- code to calculate summer mean and sd of N, P and % algal cover
trophic_group_proportions.R <-- code to calculate proportion of trophic groups in streams based on visual stomach content identification
baseline_outliers.R <-- code to determine outliers in stable isotope baselines
coupling_tp_analysis.R <-- code to calculate seasonal % terrestrial energy use and trophic position of creek chub 
isotopic_niche_analysis.R <-- code to calculate estiamtes of seasonal isotope niche space and overlap for creek chub
figure_2 <-- code to generate figure 2 (combine % terrestrial energy use, trophic position, and isotopic niche plots)
fish_size_spectrum_MLE.R <-- code to calculate size spectrum slopes using MLE methods adapted from Edwards et al., 2017, Methods in Ecology and Evolution, 8:1, 57-67

Data Files: 
fis_baseline_data.csv <-- all fish and baseline data (species, weights, stomach contents, raw stable isotope values). NA indicate information not available for given sample. 
physical_stream_conditions.csv <-- wetted width, water depth and hydraulic head measurements from each stream
summer_algae_cover.csv <-- summer % algal cover estimates
summer_nutrient_data.csv <-- summer nutrient concentrations of N and P
