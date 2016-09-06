function [ fh ] = checkPlots( avg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
plab = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Density x Length', ...
    'Weighted Fiber Count', 'Weighted Fiber Density', 'Weighted Fiber Length', 'Weighted Fiber Density x Length', ...
    'Sum of Weights', 'Weights / Count', 'Weights / Density', 'Weights / Length', ...
    'Strength of Evidence', 'Earth Movers Distance', 'Jeffery''s Divergence', 'Kullback-Leibler'};

% plot figure
fh = figure('Position', [325 25 1550 1175]);
for kk = 1:size(avg, 2)
    subplot(4, 4, kk);
    colormap('hot');
    imagesc(avg{kk}.emat.mean);
    axis('square'); axis('equal'); axis('tight');
    title([plab{kk} ' Mean']);
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    y = colorbar;
    ylabel(y, 'Strength of Connection');
end

end

