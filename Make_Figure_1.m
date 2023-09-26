clear all, warning off

% Set parameter values for the model

global D A

% Dispersal component
D = [1    3].*1e-4;

% Spatiotemporal discretisation details
T = 5000; X = 250; x = linspace(0,1,X);

% Competition component
a_ij = 0.99;
A = [1.0     a_ij; ... % Interaction matrix
    1./a_ij 1.0   ];

% Call the RD subroutine
u_example = ReactionDiffusion_two_spp(T,X);

% Plot the two species abundances at equilibrium for a particular set of parameters
figure(1), clf

subplot('position',[0.1 0.8 0.77 0.18])
hold on, box on; FS = 16;

CM = viridis(6);
plot(x,u_example(end,:,1),'linewidth',3,'color',CM(5,:))
plot(x,u_example(end,:,2),'linewidth',3,'color',CM(2,:))
xlim([-0.05 1.05]), ylim([0 1.1])
text(-0.135,1.15,'(A)','fontsize',FS)
xlabel('Location','fontsize',FS)
ylabel('Relative abundance','fontsize',FS)
L = legend('Species 1','Species 2');
set(L,'fontsize',FS,'box','off','position',[0.43 0.85 0.1 0.07])

subplot('position',[0.1 0.15 0.77 0.5])
hold on, box on

clear all, warning off

% Set parameter values for the model
global D A

% Go through a full range of parameter values to search for coexistence
DoAnalyses = 0;
if DoAnalyses == 1
    G = 100;
    T = 10000; 
    X = 30; 
    
    % Dispersal values
%     DV = [linspace(log(0.0001),log(0.05),G)' log(0.0001).*ones(G,1)];
    DV = [linspace(log(0.0001),log(0.0005),G)' log(0.0001).*ones(G,1)];

    DV = [DV(end:-1:1,[2 1]); [DV(1,1) DV(1,1)]; DV];
    eDV = exp(DV);

    % Competition strength
%     AV = linspace(log(3),log(1),G)';
    AV = linspace(log(1.15),log(1),G)';
    eAV = [[exp(AV(1:end-1)) 1./exp(AV(1:end-1))]; [1./exp(AV(end:-1:1)) exp(AV(end:-1:1))]]

    for d = 1:length(eDV)
        disp(eDV(d,:))
        for a = 1:length(eAV)
            
            D = eDV(d,:);

            A = [1.0     eAV(a,1); ...
                 eAV(a,2) 1.0   ];
            
            % Solve for equilibrium
            u_s = ReactionDiffusion_two_spp(T,X);

            % What's the coexistence situation? Zero or one indicate dominance.
            Coexist(d,a) = sum(u_s(end,:,1))./sum(sum(u_s(end,:,:)));

            keyboard
        end
    end
    
    % We assert that low abundances are equivalent to exclusion (i.e., not coexistence).
    Coexist_full = Coexist;
    Coexist(Coexist<0.01) = nan;
    Coexist(Coexist>0.99) = nan;
    
    % Coexistence is symmetrical if you define it as the minimum abundance of the two species
    Coexist = min(Coexist,1-Coexist);
%     save Coexist_1 Coexist* *AV *DV
else
    load Coexist_1
end

Coexist_full(Coexist_full<0.01) = nan;
Coexist_full(Coexist_full>0.99) = nan;

% Plot the results
FS = 16;

Coexist_full(Coexist_full<0.01) = nan;
Coexist_full(Coexist_full>0.99) = nan;
YY = eDV(:,2)./eDV(:,1);
XX = eAV(:,1)./eAV(:,2);
PP = pcolor(Coexist_full); axis tight; set(PP,'edgecolor','none')
plot([1 199],[100 100],'k--')
plot([100 100],[1 200],'k--')

CCMM = flipud(viridis(50).^0.7);
CCMM = CCMM(10:end,:);
colormap(CCMM)
cbr = colorbar;
cbr.Position = [0.2, 0.025, 0.6, 0.02]; % [x, y, width, height]
cbr.Orientation = 'horizontal';

caxis([0 1])

% xlim([90 110])
% ylim([75 125])

set(cbr,'ytick',[0 0.5 1],'yticklabel',{'0\%','50\%','100\%'})
cbr.TickLabelInterpreter = 'latex';

set(gca,'layer','top','xtick',[20:40:180],'xticklabel',[],'fontsize',18)
set(gca,'ytick',[20:40:180],'yticklabel',[])

text(-10,150,{'Species 1','superior disperser'},'rotation',90,'horizontalalignment','center','fontsize',FS)
text(-10,50,{'Species 2','superior disperser'},'rotation',90,'horizontalalignment','center','fontsize',FS)

text(50,215,{'Species 1','dominant competitor'},'horizontalalignment','center','fontsize',FS)
text(150,215,{'Species 2','dominant competitor'},'horizontalalignment','center','fontsize',FS)

text(215,190,{'$D_1 = 5$','$D_2 = 1$'},'horizontalalignment','center','fontsize',FS)
text(215,10,{'$D_1 = 1$','$D_2 = 5$'},'horizontalalignment','center','fontsize',FS)

text(215,100,{'$D_1 = 1$','$D_2 = 1$'},'horizontalalignment','center','fontsize',FS)

text(0,-10,{'$a_{12} = 0.87$','$a_{21} = 1.15$'},'horizontalalignment','left','fontsize',FS)
text(200,-10,{'$a_{12} = 1.15$','$a_{21} = 0.87$'},'horizontalalignment','right','fontsize',FS)
text(100,-10,{'$a_{ij} = 1$'},'horizontalalignment','center','fontsize',FS)

text(100,-30,{'Relative abundance of species 1'},'horizontalalignment','center','fontsize',FS+2)

text(-10,215,'(B)','horizontalalignment','center','fontsize',FS)


% Dimensions = [0 0 28 40]*0.8;
% FileName = 'Figure_1_All';
% Resolution = '-r200';
% set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
% set(gcf, 'renderer', 'painters')
% print('-dtiff',Resolution,[FileName '.tiff'])


function u = ReactionDiffusion_two_spp(T,X)

% Define the spatial and temporal domain of the PDE solution
x = linspace(0,1,X);
t = linspace(0,T,3);

% Solve the equations using matlab's PDEPE function
u = pdepe(0,@RD_eqn,@RD_ICs,@RD_BCs,x,t);

% Have we found coexistence (to within 10^6)

end

%% Define the function at this point in space/time
function [c,f,s] = RD_eqn(x,t,u,DuDx)

%% The general form of the PE equation is:
% c(x,t,u,Du/Dx) * Du/Dt = x^(-m) * D(x^m * f(x,t,u,Du/Dx))/Dx + s(x,t,u,Du/Dx)

%% Our specific form is:
% Du1/Dt = E1 D/Dx(Du1/Dx) + r*u1*(1-u1-a21*u2)
% Du2/Dt = E2 D/Dx(Du2/Dx) + r*u2*(1-u2-a12*u1)

% Try to find a combination of E and a_ij that allows coexistence
% (a12 ~= a21)

global D A
r = 1;
E = D;
a11 = A(1,1);
a22 = A(2,2);
a21 = A(2,1);
a12 = A(1,2);

c = [1; 1];
f = [E(1)*DuDx(1); E(2)*DuDx(2)];
s = [r*u(1)*(1 - a11*u(1) - a21*u(2)); ...
    r*u(2)*(1 - a22*u(2) - a12*u(1))];
end

%% Initial conditions
function u0 = RD_ICs(x)
u0 = [0.2; 0.2];
end

%% Dirichlet boundary conditions
function [pl,ql,pr,qr] = RD_BCs(xl,ul,xr,ur,t)

pl = ul-0;
ql = [0; 0];
pr = ur-0;
qr = [0; 0];

end







