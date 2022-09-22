%% Competition RDE model in 2D

function Coexist_2D_search()

for Is = 3
   
    % Load the outlines of the islands in the Kimberley archipelago
    load Kimberley_island_outlines
    
    % How many species do we want to simulate?
    NumSpp = 4;
    
    % Growth rate of all the species
    r = 1;
    
    if Is == 1
        
        % Load the appropriate island outline
        IslandOutline = StandardIslands{61};
        
        % Species interaction matrix
        A = [ 1.0000    1.0437    1.0656    1.0912
            0.9563    1.0000    1.0487    1.0588
            0.9344    0.9513    1.0000    1.0859
            0.9088    0.9412    0.9141    1.0000];
        
        % Species dispersal matrices
        D = [1.0760    2.3620    4.9864   17.7934];
        
    elseif Is == 2
        IslandOutline = StandardIslands{89};
        A = [ 1.0000    1.0738    1.0743    1.1053
            0.9262    1.0000    1.0931    1.1025
            0.9257    0.9069    1.0000    1.1211
            0.8947    0.8975    0.8789    1.0000];
        D = [0.0002    0.0129    0.0859    0.5156];
    elseif Is == 3
        IslandOutline = SS{822};
        A = [1.0000    1.0602    1.0800    1.1415
            0.9398    1.0000    1.1500    1.0750
            0.9236    0.9600    1.0000    1.1600
            0.8585    0.9250    0.8600    1.0000];
        D = [0.1855    0.2733    0.3319    1.1324];
    end
    
    %% ===== FD solution scheme =====
    
    Coexist_threshold = 1;
	a = 0; b = round(w*25);
    h = 1;
    T = 100;
    delt = 1;
    J = round((b-a)/h);
    dimJ = J+1;
    n = (dimJ)^2;
    N = round(T/delt);
    
    % Initialization
    u = zeros(n,NumSpp);
    F = zeros(n,NumSpp);
    U_grid = zeros(dimJ,dimJ,NumSpp);
    L = sparse(n,n);
    
    % Set up the domain
    [XX,YY] = meshgrid(1:dimJ,1:dimJ);
    xx = XX(:); yy = YY(:);
    
    % Fit this outline to the domain being modelled
    IslandOutline(:,1) = IslandOutline(:,1) - min(IslandOutline(:,1));
    IslandOutline(:,1) = b.*IslandOutline(:,1)./max(IslandOutline(:,1));
    IslandOutline(:,2) = IslandOutline(:,2) - min(IslandOutline(:,2));
    IslandOutline(:,2) = b.*IslandOutline(:,2)./max(IslandOutline(:,2));
    
    f_ed = find(inpolygon(xx,yy,IslandOutline(:,1),IslandOutline(:,2))==0);
    
    %% INITIAL CONDITIONS:
    N_i = 0.2;
    U0 = (N_i + XX.*0)';
    for s = 1:NumSpp
        u(:,s) = U0(:); % Convert to 1-D vector
    end
    
    % Construct matrix L (without 1/h^2 factor)
    L(1,1) = 3; L(1,2) = -3/2; L(J+1,J+1) = 6; L(J+1,J) = -3;
    L = L+sparse(2:J,3:J+1,-1,n,n);
    L = L+sparse(2:J,2:J,4,n,n);
    L = L+sparse(2:J,1:J-1,-1,n,n);
    L(1,J+2) = -3/2; L(J+1,2*J+2) = -3;
    L = L+sparse(2:J,J+3:2*J+1,-2,n,n);
    L(n-J,n-J) = 6; L(n-J,n-J+1) = -3;
    L(n,n) = 3; L(n,n-1) = -3/2;
    L = L+sparse(n-J+1:n-1,n-J+2:n,-1,n,n);
    L = L+sparse(n-J+1:n-1,n-J+1:n-1,4,n,n);
    L = L+sparse(n-J+1:n-1,n-J:n-2,-1,n,n);
    L(n-J,n-(2*J+1)) = -3; L(n,n-dimJ) = -3/2;
    L = L+sparse(n-J+1:n-1,n-2*J:n-(J+2),-2,n,n);
    L = L+sparse(J+2:n-dimJ,2*J+3:n,-1,n,n);
    L = L+sparse(J+2:n-dimJ,1:n-2*dimJ,-1,n,n);
    L = L+sparse(J+2:n-dimJ,J+2:n-dimJ,4,n,n);
    L = L+sparse(J+2:n-(J+2),J+3:n-dimJ,-1,n,n);
    L = L+sparse(J+2:dimJ:n-(2*J+1),J+3:dimJ:n-2*J,-1,n,n);
    L = L+sparse(2*J+2:dimJ:n-2*dimJ,2*J+3:dimJ:n-(2*J+1),1,n,n);
    L = L+sparse(J+3:n-dimJ,J+2:n-(J+2),-1,n,n);
    L = L+sparse(2*J+2:dimJ:n-dimJ,2*J+1:dimJ:n-(J+2),-1,n,n);
    L = L+sparse(2*J+3:dimJ:n-(2*J+1),2*J+2:dimJ:n-2*dimJ,1,n,n);
    
    % Construct matrices B1 & B2 - these represent diffusion
    for s = 1:NumSpp
        eval(['B' num2str(s) ' = sparse(1:n,1:n,1,n,n)+D(s)*L;'])
    end
    
    % Time-stepping procedure
    Diff = inf; count = 1;
    Sp_log(count,:) = sum(u);
    
    while Diff > 0.05
        
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
            %             if flagu~=0 flagu;relresu;iteru;error('GMRES did not converge'),end
        end
        
        % Apply the boundary conditions
        u(f_ed,:) = 0;
        
        count = count + 1;
        Sp_log(count,:) = sum(u);
        
        % Break the time loop if we've converged on equilibrium
        Diff = sum(abs(Sp_log(end,:)-Sp_log(end-1,:)));
        
        if min(sum(u)) < Coexist_threshold
            break
        end
        
        if count > 2500
            break
        end
    end
    
    % Did every species dominate somewhere?
    [MaxA,DomS] = max(u,[],2);
    DomS(f_ed) = [];
    H_Dom = histc(DomS,[1:NumSpp])';
    
    % Did we find coexistence?
    if Is == 1
        save Plotting_Island_1
    elseif Is == 2
        save Plotting_Island_2
    elseif Is == 3
        save Plotting_Island_3
    end
    
end





