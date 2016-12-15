function [fh] = plotMatrixOnly(inMat, crng)

% create color scale
mp = (crng(2) + crng(1)) / 2;
apts = [crng(1) mp crng(2)];
albs = {crng(1), mp, crng(2)};

% do the plot thing
colormap('hot');
imagesc(log(inMat));
axis('square'); axis('equal'); axis('tight');
caxis(crng);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);

end