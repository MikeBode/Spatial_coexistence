function Measure_coexistence_potential_N_spp()

% This is how many trials we would like to run for each island
NumRuns = 860;

load Kimberley_island_outlines_standardised StandardIslands

a = 0; b = 50; r = 1;
NumSpp = 4;
SETUP % Set up the numerical solution

% Load up the existing saved coexistence results
load Standardised_results_4_Spp Coex
COUNT = length(Coex);
% COUNT = 0; Coex = [];

% Go through all 100 islands
for IslandsTested = 1:100
    
    % Extract an island
    IslandOutline = StandardIslands{IslandsTested};
    
    % What is its area?
    PA = polyarea(IslandOutline(:,1),IslandOutline(:,2));
    f_ed = find(inpolygon(xx,yy,IslandOutline(:,1),IslandOutline(:,2))==0);
    
    Attempts = sum(Coex(:,1)==IslandsTested);
    while Attempts < NumRuns
        Attempts = sum(Coex(:,1)==IslandsTested);
        
        disp([IslandsTested Attempts])
        
        %% Construct 4 random species - using exponential distributions
        D = 2*rand; A = eye(NumSpp);
        Interactions = sort(exp(CompHeight*rand([1,NumSpp-1])),'ascend');
        A(1,:) = [1 1./Interactions];
        A(:,1) = [1 Interactions]';
        
        for s = 2:NumSpp
            
            % Dispersal ability drops geometrically
            D(s) = D(s-1)*(rand/2);
            
            % Competitive ability declines exponentially
            Interactions = sort(exp(CompHeight*rand([1,NumSpp-s])),'ascend');
            A(s,s+1:NumSpp) = 1./Interactions;
            A(s+1:NumSpp,s) = Interactions';
            
        end
        
        % Construct matrices B1 & B2 - these represent diffusion
        for s = 1:NumSpp
            eval(['B' num2str(s) ' = sparse(1:n,1:n,1,n,n) + D(' num2str(s) ')*L;'])
        end
        
        %% =-=-=-=-=- Make sure each species can persist in isolation =-=-=-=-=-
        for s = 1:NumSpp
            
            % Initial conditions for each species in isolation
            U0 = (0.25 + XX.*0)';
            u = 0.25*ones(length(U0(:)),1);
            
            eval(['B_ss = B' num2str(s) ';'])
            
            [u,Sp_log] = sub_equil_SingleSpecies(u,L,r,D,1,delt,f_ed,B_ss);
            
            if Sp_log(end) < 5e1;
                break % Don't continue if this species can't persist by itself
            end
        end
        
        %% =-=-=-=-=- Competition =-=-=-=-=-
        % Initial conditions for each species in isolation
        U0 = (1./NumSpp + XX.*0)';
        for s = 1:NumSpp
            u(:,s) = 1.*U0(:); % Convert to 1-D vector
        end
        
        [u,Sp_log] = sub_equil(u,L,r,D,A,delt,f_ed,B1,B2,B3,B4);
        
        % Store the outcomes in two variables called COEX and Coex
        COUNT = COUNT + 1;
        Coex = [Coex; IslandsTested Sp_log(end,:)];
        round(Sp_log(end,:))
        
        clearvars -except Coex Attempts CompHeight f_ed n L XX r delt COUNT IslandsTested StandardIslands xx yy Num*
    end
    
    disp(IslandsTested)
    save Standardised_results_4_Spp Coex StandardIslands
end
end





function [u,Sp_log] = sub_equil(u,L,r,D,A,delt,f_ed,B1,B2,B3,B4)

NumSpp = size(u,2);

% Time-stepping procedure
Diff = inf; count = 1;
Sp_log(count,:) = sum(u);
CoexistFlag = 1;
GMRES_flag = 0;
while Diff > 5e-2
    
    % Reaction component
    for s1 = 1:NumSpp
        F(:,s1) = u(:,s1).*r; % Population growth
        for s2 = 1:NumSpp
            F(:,s1) = F(:,s1) - A(s1,s2).*u(:,s1).*u(:,s2); % Competition
        end
    end
    u = u + delt*F;
    
    % Apply the boundary conditions
    u(f_ed,:) = 0;
    
    % Solve diffusion using GMRES
    for s = 1:NumSpp
        eval(['B = B' num2str(s) ';'])
        [u(:,s),flagu,relresu,iteru]=gmres(B,u(:,s),10,1e-8,[],[],[],u(:,s));
        if flagu~=0 flagu;relresu;iteru;disp('GMRES did not converge'), GMRES_flag = 1; break; end
    end
    
    count = count + 1;
    Sp_log(count,:) = sum(u);
    
    Diff = sum(abs(Sp_log(end,:)-Sp_log(end-1,:)));
    
    if max(u(:)) < 5e-2
        CoexistFlag = nan;
    end
    
    if GMRES_flag == 1
        CoexistFlag = nan;
        break
    end
end

end




function [u,Sp_log] = sub_equil_SingleSpecies(u,L,r,D,A,delt,f_ed,B)

NumSpp = size(u,2);

% Time-stepping procedure
Diff = inf; count = 1;
Sp_log(count,:) = sum(u);
CoexistFlag = 1;
GMRES_flag = 0;
while Diff > 5e-2
    
    % Reaction component
    for s1 = 1:NumSpp
        F(:,s1) = u(:,s1).*r; % Population growth
        for s2 = 1:NumSpp
            F(:,s1) = F(:,s1) - A(s1,s2).*u(:,s1).*u(:,s2); % Competition
        end
    end
    u = u + delt*F;
    
    % Apply the boundary conditions
    u(f_ed,:) = 0;
    
    % Solve diffusion using GMRES
    [u,flagu,relresu,iteru]=gmres(B,u,10,1e-8,[],[],[],u);
    if flagu~=0 flagu;relresu;iteru;disp('GMRES did not converge'), GMRES_flag = 1; break; end
    
    count = count + 1;
    Sp_log(count,:) = sum(u);
    
    Diff = sum(abs(Sp_log(end,:)-Sp_log(end-1,:)));
    
    if max(u(:)) < 5e-2
        CoexistFlag = nan;
    end
    
    if GMRES_flag == 1
        CoexistFlag = nan;
        break
    end
end

end








