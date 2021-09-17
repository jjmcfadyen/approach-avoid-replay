function cluster_timefrequency(filename)
% filename = .mat file containing 'subject' and 'directories' structure

addpath('/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/toolboxes/fieldtrip-20191119/');
ft_defaults;

disp(['Loading ' filename])
load(filename);
compute_timefrequency(subject,directories,waveletwidth);

end