%% dummy check of all results - not going well...
% Brent McPherson
% 20160626
% quickly load and inspect repeated data sets across subjects
%

% Franco's beatiful fuckin snowflake file structure wasted another week or
% so of compute time

% Build off aparc+aseg directly? going to do anyway...

% or make sure file loaded are correct. again.

% fine
z = feMergeRepeats('hcp', '105115', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('hcp', '105115', 'prob', '10');
fh = checkPlots(z);

% very few relative to others
z = feMergeRepeats('hcp', '110411', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('hcp', '110411', 'prob', '10');
fh = checkPlots(z);

% very few relative to others
z = feMergeRepeats('hcp', '111312', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('hcp', '111312', 'prob', '10');
fh = checkPlots(z);

% fine
z = feMergeRepeats('hcp', '113619', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('hcp', '113619', 'prob', '10');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'FP', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'FP', 'prob', '10');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'HT', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'HT', 'prob', '10');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'KK', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'KK', 'prob', '10');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'MP', 'detr', '2');
fh = checkPlots(z);

% fine
z = feMergeRepeats('stn', 'MP', 'prob', '10');
fh = checkPlots(z);

