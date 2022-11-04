function [] = DMOS_BLOOM_plot_data(day,Ehux,nano,pico,Syn,bac,EhV,...
              cil,thr,IEC,NO3,Chl,DMS,DMSPd_raw,DMSPd,DMSPp,fignum)

num_bag = size(Ehux,1);
ndays = ceil(max(day(:)));

% Color associated to each bag, following the code of Vincent et al., 2021
% color = [1 0 0; 1 0.5 0; 1 1 0; 0.8 0 0.9; 0 0 1; 0 1 1; 0 0.7 0; 0 0 0];
% faceColor = [1 1 1; 1 1 1; 1 1 1; 0.8 0 0.9; 1 1 1; 1 1 1; 0 0.7 0; 0 0 0];

% Set of colors more friendly to color-blind people
color = [0.5 0 0; 0.7 0.7 0; 0.7 0 0.7; 1 0 0; 0 0 0.5; 0 0.7 0.7; 0 0 1; 0 0 0];
faceColor = [1 1 1; 1 1 1; 1 1 1; 1 0 0; 1 1 1; 1 1 1; 0 0 1; 0 0 0];

lineType = {'o--';'s--';'d--';'o- ';'v--';'^--';'s- ';'d- '};

%..........................................................................

figure(fignum(1))
subplot(2,2,1)
ylab = 'E. huxleyi concentration [uM-N]';
DMOS_BLOOM_plot_panel(Ehux,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,2)
ylab = 'Other nanophytoplankton concentration [uM-N]';
DMOS_BLOOM_plot_panel(nano,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,3)
ylab = 'Picoeukayotic phytoplankton concentration [uM-N]';
DMOS_BLOOM_plot_panel(pico,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,4)
ylab = 'Synechococcus phytoplankton concentration [uM-N]';
DMOS_BLOOM_plot_panel(Syn,day,color,faceColor,lineType,ndays,ylab)

%..........................................................................

figure(fignum(2))
subplot(2,2,1)
ylab = 'Bacteria concentration [uM-N]';
DMOS_BLOOM_plot_panel(bac,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,2)
ylab = 'E. huxleyi virus concentration [mcp copies / L]';
DMOS_BLOOM_plot_panel(EhV,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,3)
ylab = 'Ciliates concentration [uM-N]';
DMOS_BLOOM_plot_panel(cil,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,4)
ylab = 'Thraustochytrids concentration [uM-N]';
DMOS_BLOOM_plot_panel(thr,day,color,faceColor,lineType,ndays,ylab)

%..........................................................................

figure(fignum(3))
subplot(2,2,1)
ylab = 'Nitrate concentration [uM]';
DMOS_BLOOM_plot_panel(NO3,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,2)
ylab = 'Chlorophyll a concentration [uM]';
DMOS_BLOOM_plot_panel(Chl,day,color,faceColor,lineType,ndays,ylab)

%..........................................................................

% Average values during second bloom (from day 13 on)
stp1 = 13;
EhVtot   = nanmean(EhV(1:7,stp1:end),2);
DMStot   = nanmean(DMS(1:7,stp1:end),2);
DMSPdtot = nanmean(DMSPd(1:7,stp1:end),2);
DMSPptot = nanmean(DMSPp(1:7,stp1:end),2);
Ehuxtot  = nanmean(Ehux(1:7,stp1:end),2);
Chltot   = nanmean(Chl(1:7,stp1:end),2);
DMS_DMSPp_Rat = DMStot ./ DMSPptot;

figure(fignum(4))
subplot(2,3,1)
hold on
for b=1:num_bag-1
    plot(EhVtot(b),Ehuxtot(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Average viral load [mcp copies/mL]')
ylabel('Average observed E. huxleyi concentration [uM-N]')
ylim([0 inf]);
grid on
subplot(2,3,2)
hold on
for b=1:num_bag-1
    plot(EhVtot(b),DMStot(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Avverage viral load [mcp copies/mL]')
ylabel('Average observed DMS concentration [uM-N]')
ylim([0 inf]);
grid on
subplot(2,3,3)
hold on
for b=1:num_bag-1
    plot(EhVtot(b),DMSPdtot(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Avverage viral load [mcp copies/mL]')
ylabel('Average observed DMSPd concentration [uM-N]')
ylim([0 inf]);
grid on
subplot(2,3,4)
hold on
for b=1:num_bag-1 
    plot(EhVtot(b),DMSPptot(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Average viral load [mcp copies/mL]')
ylabel('Average observed DMSPp concentration [uM-N]')
ylim([0 inf]);
grid on
subplot(2,3,5)
hold on
for b=1:num_bag-1
    plot(EhVtot(b),Chltot(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Average viral load [mcp copies/mL]')
ylabel('Average observed chlorophyll concentration [ug L{-1}]')
ylim([0 inf]);
grid on
subplot(2,3,6)
hold on
for b=1:num_bag-1 
    plot(EhVtot(b),DMS_DMSPp_Rat(b),'.','MarkerSize',30,'Color',color(b,:))
end
xlabel('Average viral load [mcp copies/mL]')
ylabel('DMS to DMSPp ratio [-]')
ylim([0 inf]);
grid on

%..........................................................................

figure(fignum(5))
hold on
for b=1:num_bag
       plot(DMSPp(b,:)+DMSPd(b,:),DMSPd(b,:),'o','MarkerSize',3,'Color',color(b,:),'MarkerFaceColor',color(b,:))
end
plot(0:800,0.2*(0:800),'r--')
plot(0:800,60*ones(1,801),'b--')

%..........................................................................

figure(fignum(6))
subplot(2,2,1)
ylab = 'Observed particulate DMSP concentration [nM]';
DMOS_BLOOM_plot_panel(DMSPp,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,2)
ylab = 'Observed (raw) dissolved DMSP concentration [nM]';
DMOS_BLOOM_plot_panel(DMSPd_raw,day,color,faceColor,lineType,ndays,ylab)
%legend({['Bag 1';'Bag 2';'Bag 3';'Bag 4';'Bag 5';'Bag 6';'Bag 7';'Fjord']},'Fontsize',15,'Location','NorthWest')
legend({'Bag 1','Bag 2','Bag 3','Bag 4','Bag 5','Bag 6','Bag 7','Fjord'},'Location','NorthWest')
set(gca,'fontsize',15); 
subplot(2,2,3)
ylab = 'Observed DMS concentration [nM]';
DMOS_BLOOM_plot_panel(DMS,day,color,faceColor,lineType,ndays,ylab)
subplot(2,2,4)
ylab = 'Observed dissolved DMSP concentration [nM]';
DMOS_BLOOM_plot_panel(DMSPd,day,color,faceColor,lineType,ndays,ylab)

%..........................................................................

figure(fignum(7))
%--------------------------------------------------------------------
subplot(3,3,1)
ylab = 'Nitrate concentration [uM]';
DMOS_BLOOM_plot_panel(NO3,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,2)
ylab = 'Chlorophyll a concentration [ug / L]';
DMOS_BLOOM_plot_panel(Chl,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,3)
ylab = 'picophytoplankton concentration [uM-N]';
DMOS_BLOOM_plot_panel(pico,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,4)
ylab = 'nanophytoplankton concentration [uM-N]';
DMOS_BLOOM_plot_panel(nano,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,5)
ylab = 'E. huxleyi concentration [uM-N]';
DMOS_BLOOM_plot_panel(Ehux,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,6)
ylab = 'Bacteria concentration [uM-N]';
DMOS_BLOOM_plot_panel(bac,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,7)
ylab = 'Emiliania huxleyi virus concentration [mcp copies / L]';
DMOS_BLOOM_plot_panel(EhV,day,color,faceColor,lineType,ndays,ylab)
%--------------------------------------------------------------------
subplot(3,3,8)
ylab = 'Infected Emiliania huxleyi cells [%]';
DMOS_BLOOM_plot_panel(IEC,day,color,faceColor,lineType,ndays,ylab)

end

