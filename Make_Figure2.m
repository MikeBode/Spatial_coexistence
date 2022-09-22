clear all

% Initialise the figure
figure(1), clf, set(gcf,'color','w'); colormap(parula); FS = 18;

% Go through the islands one-by-one
for island = 1:3
    
    % Load the dominance patterns at equilibria 
    % (These were created by 'Simulate_Island_equilibria_Fig_2.m')
    eval(['load ../Plotting_Island_' num2str(island)])
    
    % Create a subplot for each
    subplot(3,1,island), hold on
    
    % Call the plotting routine for the competition/colonisation trade-off plots
    sub_com_col_trade(A,D)
    axis square
    clearvars -except island FS
    
    if island == 1
        xlim([-0.5 3.5])
        ylim([0.92 1.08])
        text(-2.2,1.08,'(A)','fontsize',FS)
    elseif island == 2
        xlim([-6 3])
        ylim([0.875 1.1])
        text(-10,1.1,'(B)','fontsize',FS)
    elseif island == 3
        xlim([-0.75 1.5])
        ylim([0.85 1.15])
        text(-1.75,1.15,'(C)','fontsize',FS)
    end
    
    xl = xlim;
    yl = ylim;
    r = [0.1 0.9];
    set(gca,'xtick',[xl(1)+(xl(2)-xl(1))*r(1) xl(1)+(xl(2)-xl(1))*r(2)],...
        'xticklabel',{'Inferior','Superior'},...
        'ytick',[yl(1)+(yl(2)-yl(1))*r(1) yl(1)+(yl(2)-yl(1))*r(2)],...
        'yticklabel',{'Inferior','Superior'})
    
    
end

% Plot an image of each of the dominance patterns
for island = 1:3
    
    if island == 1
        w = 0.25; XX = 0.65;
        axes('position',[XX 0.65 w w]), hold on,
        axis square off vis3d; view(-30,90)
    elseif island == 2
        w = 0.25; XX = 0.65;
        axes('position',[XX 0.35 w w]), hold on, axis square off
        axis square off vis3d; view(150,90)
    elseif island == 3
        w = 0.3
        axes('position',[0.65 0.05 w w]), hold on, axis square off
        axis square off vis3d; view(110,90)
    end
    
    eval(['load ../Plotting_Island_' num2str(island)])
    
    % Apply the boundary conditions
    u(f_ed,:) = 0;
    
    % Re-order 1-D solution vectors into 2-D solution grids
    for s = 1:NumSpp
        Sp_grid(:,:,s) = reshape(u(:,s),dimJ,dimJ)';
    end
    Sp_grid(Sp_grid < 1e-2) = nan;
    
    % Plot solutions
    [YY,Dominant] = max(Sp_grid,[],3); Dominant(sum(isnan(Sp_grid),3)==NumSpp) = nan;
    
    plot(IslandOutline(:,1)+0.5,IslandOutline(:,2)+0.5,'color',0.25.*ones(1,3),'linewidth',2)
    pp = pcolor(Dominant');
    set(pp,'facealpha',0.8)
    
    CMM = viridis(4);
    colormap(CMM);
    shading flat;
    xlim([-5 dimJ+5]); ylim([-5 dimJ+5])
    
    clearvars -except island
    
end

Dimensions = [0 0 15 30];
FileName = 'Figure_2';
Resolution = '-r200';
set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print('-dtiff',Resolution,[FileName '.tiff'])
print('-depsc',Resolution,[FileName '.eps'])


function sub_com_col_trade(A,D)
cla, hold on, FS = 16; CMM = viridis(4);
for s = 1:4
    plot(log(D(s)),mean([2-A(s,:) A(:,s)']),'x','markersize',11,'linewidth',4,'color',CMM(s,:))
end
xlabel('Dispersal','fontsize',FS)
ylabel('Competition','fontsize',FS)

end

















