function plot_figures(figure_num, funcs, plot_vars, x_lim, y_lim, sub_vars, sub_vals, contour_levels, colors)
num_plots = size(sub_vals, 2);
num_funcs = size(funcs, 2);
num_rows = floor(num_plots/2);
if num_plots == 1
    num_cols = 1;
else
    num_cols = 2;
end

if num_funcs ~= size(colors, 2)
    % Adapt this for more than 7 functions
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k'];
end

figure(figure_num);
for i = 1:num_plots
    subplot(num_rows, num_cols, i);
    hold on;
    for j = 1:num_funcs
        contourSpotless(funcs(j), plot_vars(1), plot_vars(2), ...
            x_lim, y_lim, sub_vars, sub_vals(:, i), contour_levels(j), {colors(j)});
    end
end
end