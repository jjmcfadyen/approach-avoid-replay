function cluster_timefrequency(filename)
% filename = .mat file containing 'subject' and 'directories' structure

addpath('/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/toolboxes/fieldtrip-20191119/');
ft_defaults;

load(filename);
compute_timefrequency(subject,optimised_time,directories,waveletwidth);

end