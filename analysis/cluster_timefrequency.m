function cluster_timefrequency(filename)
% filename = .mat file containing 'subject' and 'directories' structure

load(filename);
compute_timefrequency(subject,optimised_time,directories,false);

end