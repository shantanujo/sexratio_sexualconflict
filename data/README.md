## This folder contains raw data files for single species trials, mixed species trials and count data.

### Explanation of column names in data files:
"`single_species_trials.csv`":\
Cage.no = Cage ID, 	sp = Species: EXS- exsulans, and TRA- traviatum, 	date = Date of trial start, 	no_male = Treatment: No. of males, 	eggs = Number of eggs, 	larvae = Number of larva, 	tand = If Tandem/Copula was observed, 	ferti = Fertility, 	female_surv = Female dead or alive, 	mat_th = Mating attempts on thorax.

"`mixed_species_trials.csv`":\
Cage.no = Cage ID ,	treatment = Treatment: E=exsulans, T= traviatum ,	date = Date of trial start ,	female_d = Female dead or alive; alive=0 ,	Tandem/Copula = If tandem was observed ,	eggs = Number of eggs ,	larvae = Number of larva ,	female_sp = Female species ,	con_male = Number of conspecific males ,	hetero_male = Number of heterospecific males ,	date_j = Day of year, ferti = Fertility ,	th_con = Mating attempts by conspecific males ,	th_het = Mating attempts from heterospecific males.

"`sex_ratio_counts.csv`":\
date = Date of survey, 	sp = Species: EXS= exsulans, and TRA= traviatum, 	M = No. of males, 	F = No. of females, 	tandem = No. of tandems, 	tot = Total no. of individuals := M + F + 2*tandem, 	date_j = Day of Year, 	year = Year, males = Total no. of males := M + tandem, 	females = Total no. of females := F + tandem, 	ratio = Total no. of males/ females.
