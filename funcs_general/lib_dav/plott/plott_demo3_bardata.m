%% Bar Graphs - Short demo
% Plot some random data and draw arbitrary significance stars



%% Example 1 - 1D barplot, 1d h matrix
y = randn(1,5);         
errY = 0.1.*y;          % 10% error
h = [0 1 0 0 1];
h = bar_errsig(h,errY, y,'g');% Plot with errorbars
ylim([-2 2.5]);

%% Example 2 - 2D barplot, 1d h matrix, with some optional arguments
% (i.e. for plotting anova output)
y = abs(randn(3,4));        
errY = 0.1.*y;          % 10% error
h = [0 1 0];
h = bar_errsig(h,errY, y,'mysymbol','#','star_shift_factor',1.2);% Plot with errorbars
ylim([0 3]);

%% Example 3 - 2D barplot, 2d h matrix
y = abs(randn(3,4));        
errY = 0.1.*y;          % 10% error
h = [0 1 0 1; 0 1 1 1; 1 0 1 0];
out = bar_errsig(h,errY, y);% Plot with errorbars
ylim([0 3]);