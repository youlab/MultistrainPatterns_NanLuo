%% Initialization
nt = totalt / dt;  % number of time steps
dx = L / (nx - 1); dy = dx;
x  = linspace(-L/2, L/2, nx);
y  = linspace(-L/2, L/2, ny);
[xx, yy] = meshgrid(x, y);
rr = sqrt(xx .^ 2 + yy .^ 2);

P = cell(3, 1); P(:) = {zeros(nx, ny)};   % Pattern (P = 1 inside colony; P = 0 ouside colony)
C = cell(3, 1); C(:) = {zeros(nx, ny)};   % Cell density
N = zeros(nx, ny) + N0;                   % Nutrient

for j = 1 : 3
    P{j}(rr <= r0) = 1;
    C{j}(P{j} == 1) = C0 * initialFract(j) / (sum(P{j}(:)) * dx * dy);
end
C_pre = C;

% branch density & width of mixed colonies
branchDensity = initialFract * Densities';
branchWidth   = initialFract * Widths';

ntips = ceil(2 * pi * r0 * branchDensity); % branch number
ntips = max(ntips, 2);  % initial branch number cannot be less than 2
% x/y coordinates of every tip of every species at every timestep 
Tipx = cell(3, 1); Tipx(:) = {zeros(ntips, totalt / dt_updatebranch + 2)};
Tipy = cell(3, 1); Tipy(:) = {zeros(ntips, totalt / dt_updatebranch + 2)};

CooperMot = zeros(ntips, 3); % cooperative motility contributed by each branch of each species

theta = linspace(0, 2 * pi, ntips + 1)' + rand * 0;
theta = theta(1 : ntips) * ones(1, 3);  % growth directions of each branch of each species
delta = linspace(-1, 1, 201) * pi;

[MatV1N,MatV2N,MatU1N,MatU2N] = diffusionADI(dx,dy,nx,ny,dt,DN); % for solving diffusion equation using ADI method

aCs = cell(3, 1);
hs  = cell(3, 1);
localFracts = cell(3, 1);

BiomassV = zeros(nt + 2, 3); 
BiomassV(1, :) = [sum(C{1},'all'); sum(C{2},'all'); sum(C{3},'all')];    

%% Growth iterations

for i = 0 : nt    
    
    fN = N ./ (N + 1) .* (1 - (C{1} + C{2} + C{3}));
    fN(fN < 0) = 0;
    
    % -------------------------------------
    % Cell growth
    for j = 1 : 3  % j: index of species
        aCs{j} = aCs_act(j) * ones(nx, ny); % growth rate as a variable of time and location
        aCs{j}(N > N_upper | N < N_lower) = max(aCs_act); % if nutrient is not within the range grow at full speed
        C{j} = C{j} + aCs{j} .* fN .* C{j} * dt;
    end
    
    % -------------------------------------
    % Nutrient distribution
    dN = - bN * fN .* (aCs{1} .* C{1} + aCs{2} .* C{2} + aCs{3} .* C{3}); % Nutrient consumption
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N); N = (MatV2N * NV) / MatU2N; % Nutrient diffusion

    % -------------------------------------
    % Cell allocation and branch extension

    if mod(i, dt_updatebranch / dt) == 0
    
        cellGrowth = cell(3, 1);

        for j = find(initialFract > 0) 
            % get cell growth of each species
            cellGrowth{j} = (C{j} - C_pre{j}) * dx * dy; % total cell growth  
            
            % save the locations of branch tips
            ib = i / (dt_updatebranch / dt) + 2;
            Tipx{j}(:, ib) = Tipx{j}(:, max(1, ib - 1));
            Tipy{j}(:, ib) = Tipy{j}(:, max(1, ib - 1));         
            
            % calculate nutrient-dependent cooperative motility coefficient
            hs{j} = hs_act(j) * ones(nx, ny);
            hs{j}(N > N_upper | N < N_lower) = 0; % swarm only when nutrient is within the range
            
            % local fractions of each species
            localFracts{j} = C{j} ./ (C{1} + C{2} + C{3}); 
        end
    
        for j = find(initialFract > 0)
            
            % update branch width and density with time  
            nn = ntips;
            tipFracts = zeros(nn, 3);
            for jj = find(initialFract > 0) 
                tipFracts(:, jj) = interp2(xx, yy, localFracts{jj}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib));
                % the fraction of each strain at each branch tip
            end
            branchDensity = tipFracts * Densities';
            branchWidth   = tipFracts * Widths';
            
            % branch width cannot be larger than the inoculum
            maxbranchWidth = sqrt(Tipx{j}(:,ib) .^ 2 + Tipy{j}(:,ib) .^ 2) + r0;
            branchWidth(branchWidth > maxbranchWidth) = maxbranchWidth(branchWidth > maxbranchWidth);
    
            % get cooperative motility contributed by each branch of each species        
            temp = zeros(nx, ny);
            for jj = find(initialFract > 0) 
                temp = temp + hs{jj} .* cellGrowth{jj}; % summing up cooperative motility contributed by all species
            end
            for k = 1 : nn
                d = sqrt((Tipx{j}(k,ib) - xx) .^ 2 + (Tipy{j}(k,ib) - yy) .^ 2);
                ind = d <= branchWidth(k)/2; % the tip of the kth branch of the jth species
                CooperMot(k,jj) = sum(temp .* ind,'all'); 
                % cooperative motility contributed by all species within the kth branch tip of the jth species
            end
    
            % extension rate of each branch  
            dl = (gs(j) + rs(j) .* sum(CooperMot(1:nn,:), 2))./ branchWidth;
            if i == 0; dl = 0.5; end
            
            % make slower species follow faster species
            [Tipx{j}, Tipy{j}] = tiptracking(Tipx, Tipy, ib, dl, theta, delta, nn, xx, yy, N, j, noiseamp);
    
            % branch bifurcation
            % a branch bifurcates if there is no other branch tips within the radius of R
            R = 3/2 ./ branchDensity;  
            Tipx_new = Tipx{j}; Tipy_new = Tipy{j}; 
            theta_new = theta(:,j); dl_new = dl;
            for k = 1 : nn
                dist2othertips = sqrt((Tipx_new(:,ib) - Tipx{j}(k,ib)) .^ 2 + (Tipy_new(:,ib) - Tipy{j}(k,ib)) .^ 2);
                dist2othertips = sort(dist2othertips);
                if dist2othertips(2) > R(k)
                    nn = nn + 1;
                    Tipx_new(nn,:) = Tipx{j}(k,:);
                    Tipy_new(nn,:) = Tipy{j}(k,:);
                    for jk = 1 : 3
                        if jk ~= j && size(Tipx{jk}, 1) < nn
                            Tipx{jk}(nn,:) = Tipx{jk}(k,:);
                            Tipy{jk}(nn,:) = Tipy{jk}(k,:);
                        end
                    end
                    Tipx_new(nn,ib) = Tipx{j}(k,ib) + 0.7654*dl(k) * sin(theta(k,j) + 0.625 * pi); % splitting the old tip to two new tips
                    Tipy_new(nn,ib) = Tipy{j}(k,ib) + 0.7654*dl(k) * cos(theta(k,j) + 0.625 * pi); % numbers are to ensure extension = dl & separation angle = 90 after splitting
                    Tipx_new(k,ib) = Tipx_new(k,ib) + 0.7654*dl(k) * sin(theta(k,j) - 0.625 * pi);
                    Tipy_new(k,ib) = Tipy_new(k,ib) + 0.7654*dl(k) * cos(theta(k,j) - 0.625 * pi);
                    dl_new(nn) = dl(k) / 2;
                    dl_new(k)  = dl(k) / 2;
                    theta_new(nn) = theta(k,j);
                    branchWidth(nn) = branchWidth(k);
                end
            end
            ntips = nn;
            if nn > size(Tipx{j}, 1)
                nz = nn - size(Tipx{j}, 1);
                theta = [theta; zeros(nz, 3)];
            end
            Tipx{j} = Tipx_new; Tipy{j} = Tipy_new; 
            theta(:,j) = theta_new; dl = dl_new; 
            
            % tips cannot go out of the petri dish
            idx = abs(Tipx{j}(:, ib)) > L/2 | abs(Tipy{j}(:, ib)) > L/2;
            Tipx{j}(idx, ib) = Tipx{j}(idx, ib - 1);
            Tipy{j}(idx, ib) = Tipy{j}(idx, ib - 1);
    
            % fill the width of the branches
            for k = 1 : nn
                d = sqrt((Tipx{j}(k,ib) - xx) .^ 2 + (Tipy{j}(k,ib) - yy) .^ 2);
                P{j}(d <= branchWidth(k)/2) = 1;
            end
    
            % reallocate biomass
            remainCapacity = 1 - C_pre{1} - C_pre{2} - C_pre{3}; % remaining cell capacity
            remainCapacity(P{j} == 0) = 0;   % no capacity outside the colony
            C{j} = C_pre{j} + sum(cellGrowth{j}(:)) / sum(remainCapacity(:)) * remainCapacity / (dx * dy);
    
            % plot
            if mod(i, round(plotgap/dt)) == 0; plotting; end
            if saveintermediateimages && mod(i, 100) == 0; MakeFigure_cooperationregions; end
    
        end
    
        C_pre = C;

    end

    BiomassV(i + 2, :) = [sum(C{1},'all'); sum(C{2},'all'); sum(C{3},'all')];
    
end

Biomass = [sum(C{1},'all'); sum(C{2},'all'); sum(C{3},'all')];
