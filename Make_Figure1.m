clear all, warning off

% Set parameter values for the model
global D A
D = [0.05; 0.0003]; % Dispersal parameters
a_ij = 3;
A = [1.0     a_ij; ... % Interaction matrix
    1./a_ij 1.0   ];

% Spatiotemporal discretisation details
T = 50; X = 100; x = linspace(0,1,X);

% Call the RD subroutine
u_example = ReactionDiffusion_two_spp(T,X);

% Plot the two species abundances at equilibrium for a particular set of parameters
figure(1), clf, FS = 16;
CM = viridis(6);
subplot(3,1,1), hold on
plot(x,u_example(end,:,1),'linewidth',3,'color',CM(5,:))
plot(x,u_example(end,:,2),'linewidth',3,'color',CM(2,:))
ylim([0 0.75])
xlim([-0.1 1.1])
text(-0.1,0.95,'(A)','fontsize',FS)
xlabel('Location','fontsize',FS)
ylabel('Relative abundance','fontsize',FS)
L = legend('Species 1','Species 2');
set(L,'fontsize',FS,'box','off','location','northoutside')

% Now go through a full range of parameter values to search for coexistence
DoAnalyses = 1
if DoAnalyses == 1
    T = 100; X = 40;
    
    % Search through all these values
    DV = linspace(-7,-4,201);
    eDV = exp(DV);
    AV = linspace(0.4,1.5,200);
    eAV = exp(AV);
    
    for d = 1:length(DV)
        for a = 1:length(AV)
            
            disp([DV(d) AV(a)])
            
            D = [0.05; eDV(d)];
            a_ij = eAV(a);
            A = [1.0     a_ij; ...
                1./a_ij 1.0   ];
            
            % Solve for equilibrium
            u_s = ReactionDiffusion_two_spp(T,X);
            
            % What's the coexistence situation? Zero or one indicate dominance.
            Coexist(d,a) = sum(u_s(end,:,1))./sum(sum(u_s(end,:,:)));
        end
    end
    
    % We assert that low abundances are equivalent to exclusion (i.e., not coexistence).
    Coexist_full = Coexist;
    Coexist(Coexist<0.025) = nan;
    Coexist(Coexist>0.975) = nan;
    
    % Coexistence is symmetrical if you define it as the minimum abundance of the two species
    Coexist = min(Coexist,1-Coexist);
    save Coexist Coexist* AV DV
else
    load Coexist Coexist AV DV
end

% Plot the results
subplot(3,1,[2 3]), cla, hold on, box on, axis square
PP = pcolor(AV,DV,Coexist); axis tight; set(PP,'edgecolor','none')
colormap(flipud(viridis(50).^0.75)), colorbar

text(0.25,-4.05,'(B)','fontsize',FS)
xlabel('Competitive superiority of species 2','fontsize',FS)
ylabel('Relative dispersal ability of species 2','fontsize',FS)
title('Abundance of least common species','fontsize',FS)

XT = [0.5 1:5];
YT = [0.001 0.002 0.005 0.01];
YTL = 0.05./YT;
set(gca,'layer','top','xtick',log(XT),'xticklabel',XT), caxis([0 0.5])
set(gca,'ytick',log(YT),'yticklabel',YTL)

Dimensions = [0 0 20 30];
FileName = 'Figure_1';
Resolution = '-r200';
set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print('-dtiff',Resolution,[FileName '.tiff'])
print('-depsc',Resolution,[FileName '.eps'])


function u = ReactionDiffusion_two_spp(T,X)

% Define the spatial and temporal domain of the PDE solution
x = linspace(0,1,X);
t = linspace(0,T,10);

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







