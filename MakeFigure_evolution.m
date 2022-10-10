% making figures for Fig 6A
figure(2); clf
set(gcf, 'position', [320.3333  345.6667  185.3333   82.6667])

Frat = BiomassDynamics ./ sum(BiomassDynamics, 2); 
area(linspace(0, ndays, (nt + 1) * ndays + 1), Frat)

color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
colororder([color1; color2; color3]/255)
ylim([0 1]); xlim([0 ndays])
xlabel 'Days'
ylabel 'Fraction'
set(gca, 'TickLength', [0,0])
set(gca, 'FontSize', 8)
set(gca, 'LabelFontSizeMultiplier', 10/8)
set(gca, 'xtick', 0 : ndays)
set(gca, 'ytick', [0 1])
set(gca, 'TitleFontWeight', 'bold')
set(gca, 'XGrid', 'on')
set(gca, 'Layer', 'top')
set(gca, 'GridAlpha', 0.5)

fig = gcf; fig.PaperPositionMode = 'auto';
print([filename '.svg'],'-dsvg')