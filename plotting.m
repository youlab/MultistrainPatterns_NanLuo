
figure(1);
set(gcf,'position',[213  73  966  538]);
iz = 1 : 2 : nx; % plot less cols & rows to save time

% Plot each species
subplot(2, 3, j)
    hold off; pcolor(xx(iz, iz), yy(iz, iz), C{j}(iz, iz));
    shading interp; axis equal; caxis([-max(C{j}(:)) max(C{j}(:))])
    axis([-L/2 L/2 -L/2 L/2]); hold on
    set(gca,'YTick',[], 'XTick',[])
    plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
    title(speciesName{j})
    drawnow

% Plot all species
if j == find(initialRatio,1,'last')

color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
Ctotal = 5*C{1} + C{2} + C{3};
p1 = 5*C{1}./Ctotal; p1(isnan(p1)) = 0; 
p2 = C{2}./Ctotal; p2(isnan(p2)) = 0;
p3 = C{3}./Ctotal; p3(isnan(p3)) = 0;

subplot(2, 3, 4) % total cell density
    hold off; pcolor(xx(iz, iz), yy(iz, iz), Ctotal(iz, iz));
    shading interp; axis equal; caxis([-max(Ctotal(:)) max(Ctotal(:))])
    axis([-L/2 L/2 -L/2 L/2]); colormap('gray'); hold on
    set(gca,'YTick',[], 'XTick',[])
    plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
    title('Total cell density')
subplot(2, 3, 5) % show each species by color
    ColorMap = MarkMixing_3color(color1, color2, color3, p1, p2, p3);
    hold off; surf(xx(iz, iz), yy(iz, iz), ones(size(xx(iz, iz))), ColorMap(iz, iz, :))
    view([0, 0, 1]); shading interp; axis equal; box on
    axis([-L/2 L/2 -L/2 L/2]);
    set(gca,'YTick',[], 'XTick',[])
    title('All')
subplot(2, 3, 6) % line graph of cell densities
    yyaxis left; hold off
    mid = (nx + 1) / 2;
    plot(x(mid:end), C{1}(mid:end,mid), '-', 'color', color1/255, 'linewidth', 2); hold on
    plot(x(mid:end), C{2}(mid:end,mid), '-', 'color', color2/255, 'linewidth', 2);
    plot(x(mid:end), C{3}(mid:end,mid), '-', 'color', color3/255, 'linewidth', 2);
    plot(x(mid:end), Ctotal(mid:end,mid), 'k-', 'linewidth', 2)
    ylabel 'Cell densities';
    yyaxis right; hold off
    plot(x(mid:end), N(mid:end,mid), '-', 'color', [0.7,0.7,0.7], 'linewidth', 2); ylim([0 N0])
    ylabel 'Nutrient'
    xlabel 'Distance from center'
drawnow

end
    