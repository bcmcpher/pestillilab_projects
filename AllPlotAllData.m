
%initialize variable structure
load('/N/dc2/projects/lifebid/HCP/Dan/LifeFiberWeightOUT.mat')
load('/N/dc2/projects/lifebid/HCP/Dan/LifeFiberCountOUT.mat')
load('/N/dc2/projects/lifebid/HCP/Dan/LifeFiberLengthVol.mat')
colors=[.8 .2 .8; .8 .6 .8; .2 .8 .8; .6 .8 .8; .2 .8 .3; .6 .2 .3; .8 .8 .2; .8 .8 .6];

savedir='/N/dc2/projects/lifebid/HCP/Dan/';

curCount=zeros(1,length(AmalgumVLWeightOut)/10);

dataDim=size(FibLengthVolStruc);

dataStruc=[];

for iconditions = 1:dataDim(1)/10
    for ifibers=1:(dataDim(2)-1)
       dataStruc.countData{iconditions,ifibers}={};
       dataStruc.weightData{iconditions,ifibers}={};
       dataStruc.lengthData{iconditions,ifibers}={};
       dataStruc.volData{iconditions,ifibers}={};
    end
end
        
        
%stupid hardcode because subject 5 in Stanford is not ready yet
Hardstop=840;


dataNames=[];

for iconditions = 1:(length(AmalgumVLWeightOut)/10)
underscoreIndex=strfind(AmalgumVLWeightOut{((iconditions*10)-1),1}, '_');
    dataNames{iconditions}=AmalgumVLWeightOut{((iconditions*10)-1),1}(1:underscoreIndex(end)-1);
end

for iconditions = 1:Hardstop/10
%     underscoreIndex=strfind(AmalgumVLWeightOut{((iconditions*10)-1),1}, '_');
%     dataNames{iconditions}=AmalgumVLWeightOut{((iconditions*10)-1),1}(1:underscoreIndex(end)-1);
    
    for ifibers=1:(dataDim(2)-1)
        fiberWeightSum=sum([AmalgumVLWeightOut{((iconditions*10)-9):iconditions*10,ifibers+1}]);
        fiberCountSum=sum([AmalgumVLOut{((iconditions*10)-9):iconditions*10,ifibers+1}]);
        fiberVolumeSum=0;
        fiberLengthSum=0;
        for iten=0:9
         
        fiberVolumeSum=fiberVolumeSum+([FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.volume]);
        fiberLengthSum=fiberLengthSum+([FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.length]);
        
        %create vector to hold all 10 reconstructions for sd and sem
        %calculation
        if isempty (FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.volume);
            fiberVolumeVec(iten+1)=0;
            
        else
        fiberVolumeVec(iten+1)=FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.volume;
        end
        if isempty (FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.length);
            fiberLengthVec(iten+1)=0;
        else
        fiberLengthVec(iten+1)=FibLengthVolStruc{((iconditions*10)-iten),ifibers+1}.length;
        end
        fiberCountVec(iten+1)=AmalgumVLOut{((iconditions*10)-iten),ifibers+1};
        fiberWeightVec(iten+1)=AmalgumVLWeightOut{((iconditions*10)-iten),ifibers+1};
        end
        
        dataStruc.countData{iconditions,ifibers}.fiberCountSEM=std(fiberCountVec)/(sqrt(10));
      
            
            fiberLengthVec=fiberLengthVec(find(fiberLengthVec>0));
            fiberVolumeVec=fiberVolumeVec(find(fiberVolumeVec>0));
          if ifibers== or(or(21,22),or(7,8))
            fiberVolumeVec
            fiberLengthVec
        end
        dataStruc.weightData{iconditions,ifibers}.fiberWeightSEM=std(fiberWeightVec)/(sqrt(10));
        dataStruc.lengthData{iconditions,ifibers}.fiberLengthSEM=std(fiberLengthVec)/(sqrt(length(fiberLengthVec)));
        dataStruc.volData{iconditions,ifibers}.fiberVolumeSEM=std(fiberVolumeVec)/(sqrt(length(fiberVolumeVec)));
        
        dataStruc.countData{iconditions,ifibers}.fiberCountAVG=mean(fiberCountVec);
        dataStruc.weightData{iconditions,ifibers}.fiberWeightAVG=mean(fiberWeightVec);
        dataStruc.lengthData{iconditions,ifibers}.fiberLengthAVG=mean(fiberLengthVec);
        dataStruc.volData{iconditions,ifibers}.fiberVolumeAVG=mean(fiberVolumeVec);
        
        dataStruc.countData{iconditions,ifibers}.fiberCountSum=fiberCountSum;
        dataStruc.weightData{iconditions,ifibers}.fiberWeightSum=fiberWeightSum;
        dataStruc.lengthData{iconditions,ifibers}.fiberLengthSum=fiberLengthSum;
        dataStruc.volData{iconditions,ifibers}.fiberVolumeSum=fiberVolumeSum;
        
    end
      
    fprintf ('\n condition %i ',iconditions)
    
end

dataStrucCopy=dataStruc;

%previously
%1 = 'L_Arcuate_Posterior'
%2 = 'R_Arcuate_Posterior'
%3 = 'L_VOF'
%4 = 'R_VOF'
%5 = 'Left Arcuate'
%6 = 'Right Arcuate'
%7 = 'LeftLatVPF'
%8 = 'LeftMedVPF'
%9 = 'RightLatVPF'
%10 ='RightMedVPF'
%11 = 'Left MdLF'
%12 = 'Right MdLF'
%13 = 'Left Thalamic Radiation'
%14 = 'Right Thalamic Radiation'
%15 = 'Left Corticospinal'
%16 = 'Right Corticospinal'
%17 = 'Left Cingulum Cingulate'
%18 = 'Right Cingulum Cingulate'
%19 = 'Left Cingulum Hippocampus'
%20 = 'Right Cingulum Hippocampus'
%21 = 'Callosum Forceps Major'
%22 = 'Callosum Forceps Minor'
%23 = 'Left IFOF'
%24 = 'Right IFOF'
%25 = 'Left ILF'
%26 = 'Right ILF'
%27 = 'Left SLF'
%28 = 'Right SLF'
%29 = 'Left Uncinate'
%30 = 'Right Uncinate'

%desired order
% [ 5 6 27 28 17 18 15 16 21 22 13 14 29 30 19 20 23 24 25 26 11 12 3 4 7 9
% 8 10]

%stupid hardcoding of fiber names because I didn't come up with a standard
%naming pattern
fiberNamesVec={'Posterior Arcuate', 'VOF', 'Arcuate', 'Lateral VPF', 'Medial VPF', 'MdLF', 'Thalamic Radiation', 'Corticospinal', 'Cingulum Cingulate', 'Cingulum Hippocampus', 'Callosum Forceps Major/Minor', 'IFOF', 'ILF', 'SLF', 'Uncinate'};

fiberSequence=[ 5 6 27 28 17 18 15 16 21 22 13 14 29 30 19 20 23 24 25 26 11 12 3 4 7 9 8 10];
nameSequence=[3 14 9 8 11 7 15 10 12 13 6 2 4 5];


for Inames=1:length(nameSequence)
    namesVec1{Inames}=fiberNamesVec{nameSequence(Inames)};
end


%flag Maker
subjectHold=[0 (1:8)*12];
for isubj=1:8
    for iiter=(subjectHold(isubj)+1):(subjectHold(isubj+1))
       
    SubjectFlags(iiter)=isubj;
    end
end

for ilmax=0:5
    for trythis=( ilmax+1):6:(96)
   LmaxFlags(trythis)=ilmax;
       
    end
end
LmaxFlags=LmaxFlags+1;

for isubjects=1:8
    subjectIndexes=find(SubjectFlags==isubjects);
    subjectRunSize=length(subjectIndexes);
    probDetFlag(subjectIndexes(1):subjectIndexes(subjectRunSize/2))=1;
    probDetFlag(subjectIndexes((subjectRunSize/2)+1):subjectIndexes(end))=2;
end

    groupFlag(1:48)=1;
    groupFlag(49:96)=2;

c1 = colormap(parula(64));
c2 = colormap(autumn(64));

c = [ c2([49 40 28 4],:);c1([4 15 27 31],:)];

markerVec={'s' 'o'};

delta = 0.8/(length(unique(SubjectFlags))-1);
jitter = [-3.5, -2.5, -1.5, -0.5, +0.5, +1.5, +2.5, +3.5]*delta;

a = 0.5;


% plotting

%% Count

close all
figure
for idata=1:84
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
        semilogx(dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG - dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountSEM ; ...
            dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG + dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Count - Major human whith matter tracts');
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
    'xlim',[1 4096], 'xtick',[1 4 16 64  256 1024 4096], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])
figureNamed=(strcat(savedir,'CountPlot_AllLmax.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')



for ilmax=1:6
lmaxfives=find (LmaxFlags==ilmax);
close all
figure
for idata=lmaxfives(1:14)
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
        semilogx(dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG - dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountSEM ; ...
            dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountAVG + dataStruc.countData{idata,fiberSequence(Ifibers)}.fiberCountSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax));
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
    'xlim',[1 4096], 'xtick',[1 4 16 64  256 1024 4096], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'CountPlot_Lmax',num2str(2*ilmax),'.png'))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')
end

%% Weight
close all
figure
for idata=1:84
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
  semilogx(dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG - dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightSEM ; ...
            dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG + dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Weight - Major human whith matter tracts');
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
    'xlim',[.01 100], 'xtick',[.1 1 10 100], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'WeightPlot_AllLmax.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')



for ilmax=1:6
lmaxfives=find (LmaxFlags==ilmax);
close all
figure
for idata=lmaxfives(1:14)
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
  semilogx(dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG - dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightSEM ; ...
            dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightAVG + dataStruc.weightData{idata,fiberSequence(Ifibers)}.fiberWeightSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Weight - Major human whith matter tracts at L_max',num2str(2*ilmax));
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
    'xlim',[.01 100], 'xtick',[.1 1 10 100], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'WeightPlot_Lmax',num2str(2*ilmax),'.png'))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')
end

%% Length
close all
figure
for idata=1:84
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
        plot(dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        plot([dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG - dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthSEM ; ...
            dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG + dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
              
    end
end

fh=gcf;
fh.Name = strcat('Length - Major human whith matter tracts');
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
    'xlim',[0 180], 'xtick',[0:20:180], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'LengthPlot_AllLmax.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')



for ilmax=1:6
lmaxfives=find (LmaxFlags==ilmax);
close all
figure
for idata=lmaxfives(1:14)
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
        plot(dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        plot([dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG - dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthSEM ; ...
            dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthAVG + dataStruc.lengthData{idata,fiberSequence(Ifibers)}.fiberLengthSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Length - Major human whith matter tracts at L_max',num2str(2*ilmax));
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
     'xlim',[0 180], 'xtick',[0:20:180], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'LengthPlot_Lmax',num2str(2*ilmax),'.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')
end






%% Volume

close all
figure
for idata=1:84
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
       semilogx(dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG - dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeSEM ; ...
            dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG + dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
    end
end

fh=gcf;
fh.Name = strcat('Volume - Major human whith matter tracts');
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
     'xlim',[16 65536], 'xtick',[ 16 64 256 1024 4096 16384 65536 ], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'VolumePlot_AllLmax.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')



for ilmax=1:6
lmaxfives=find (LmaxFlags==ilmax);
close all
figure
for idata=lmaxfives(1:14)
    for Ifibers=1:length (fiberSequence)
        x = ( Ifibers+ jitter(SubjectFlags(idata)));
        % prob
       semilogx(dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG,x,markerVec{probDetFlag(idata)},'markerfacecolor',c(SubjectFlags(idata),:),  'markeredgecolor','k','linewidth',0.5,'markersize',10)
        hold on
        
        semilogx([dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG - dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeSEM ; ...
            dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeAVG + dataStruc.volData{idata,fiberSequence(Ifibers)}.fiberVolumeSEM ],  [x; x], '-','color',[a a a],'linewidth',2)
        
        
    end
end

fh=gcf;
fh.Name = strcat('Volume - Major human whith matter tracts at L_max',num2str(2*ilmax));
set(fh,'Position',[0,0,700,1500]);

%title(strcat('Count - Major human whith matter tracts at L_max',num2str(2*ilmax))')


set(gca, ...
     'xlim',[16 65536], 'xtick',[ 16 64 256 1024 4096 16384 65536 ], ...
    'ylim',[0 28], 'ytick',[1.5:2:28], ...
    'yticklabel',namesVec1, ...
    'tickdir','out', ...
    'box','off', ...
    'ticklen',[0.01 .01])

figureNamed=(strcat(savedir,'VolumePlot_Lmax',num2str(2*ilmax),'.png'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-dpng', figureNamed ,'-r100')
end

%%
for isubjectFlags=1:8
subjectIndexes(isubjectFlags)=find(isubjectFlags==SubjectFlags);
end




%% end
