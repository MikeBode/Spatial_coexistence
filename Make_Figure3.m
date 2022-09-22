clear all

% Load the standardised island outlines and coexistence results
load Standardised_results_4_Spp Coex StandardIslands

% Calculate the perimeter (the area is constant)
for i = 1:100
	ThisIslandPolygon = StandardIslands{i};

    Metric(i) = Code_Dylan_Maher(ThisIslandPolygon);
    Metric_2(i) = Code_William_Collum(ThisIslandPolygon);
end

% Set up the color scheme (PERIMETER)
for i = 1:100 
    SDC(i) = sum(Metric < Metric(i))+1; 
end
CL = flipud(viridis(100));

% Count how many times coexistence occurred on each island
MinCoexist = 50;
for i = 1:100
    ThisIsland = find(Coex(:,1) == i);
    PersistThisIsland = Coex(ThisIsland,2:end) > MinCoexist;
    H(i) = mean(sum(PersistThisIsland,2));
end

[R,I] = sort(H,'ascend');

% Plot the islands in order of coexistence
figure(1), clf; count = 0; d = 0.12; w = 0.275; hh = 0.975;
subplot('position',[0.07 0.07 0.5 0.85]), hold on, set(gcf,'color','w'); axis off equal
for i = 1:1:100
    Y = 10-ceil(i/10); Y = -Y.*1.8;
    X = mod((i-1),10); X = X.*1.1;

    ThisIsland = StandardIslands{I(i)};
    if length(ThisIsland) > 0
        ThisIsland(:,1) = ThisIsland(:,1)-min(ThisIsland(:,1));
        ThisIsland(:,1) = ThisIsland(:,1)./max(ThisIsland(:,1));
        ThisIsland(:,2) = ThisIsland(:,2)-min(ThisIsland(:,2));
        ThisIsland(:,2) = ThisIsland(:,2)./max(ThisIsland(:,2));
        pp = patch(X+ThisIsland(:,1),Y+ThisIsland(:,2),'k');
        set(pp,'facealpha',0.7,'facecolor',CL(SDC(I(i)),:),'edgecolor',CL(SDC(I(i)),:).^2,'linewidth',1.5)
        count = count + 1;
        
    end
    drawnow
end

annotation('arrow',[0.05 0.05],[0.2 0.85],'linewidth',2)
text(-1.5,-8,'Increasing coexistence','fontsize',17,'color','k','horizontalalignment','center','rotation',90)
ylim([-17 2])
xlim([-0.5 11.5])

colormap(CL); 
cb = colorbar('SouthOutside');
cb.Ticks = [0.15 0.85];
cb.TickLabels = {'Simple geometry';'Complex geometry'};
cb.FontSize = 17;
cb.TickLabelInterpreter = 'Latex';

MS = 10; FS = 16;

text(-1.4,1.6,'(A)','fontsize',FS+2)

subplot('position',[0.68 hh-d-w w w]), hold on, set(gca,'linewidth',1.5); box on, axis square
Metric = Metric./max(Metric);
plot(Metric,H,'.','markersize',MS,'color',0.3.*ones(1,3))
xlabel('Complexity metric \#1','fontsize',FS)
ylabel('Coexisting species','fontsize',FS)
axis square
box on
xlim([0.2 1])
text(0.23,2.44,'(B)','fontsize',FS+2)

subplot('position',[0.68 hh-2*d-2*w w w]), hold on, set(gca,'linewidth',1.5); box on, axis square
Metric_2 = Metric_2./max(Metric_2);
plot(Metric_2,H,'.','markersize',MS,'color',0.3.*ones(1,3))
xlabel('Complexity metric \#2','fontsize',FS)
ylabel('Coexisting species','fontsize',FS)
axis square
box on
xlim([0.2 1])
text(0.23,2.44,'(C)','fontsize',FS+2)


Dimensions = [0 0 30 30];
FileName = 'Figure_3';
Resolution = '-r200';
set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print('-dtiff',Resolution,[FileName '.tiff'])
print('-depsc',Resolution,[FileName '.eps'])


