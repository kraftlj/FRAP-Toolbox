%%% FRAP Toolbox
%%
%     FRAP Toolbox is designed to be a modular software program designed
%     for the purposes of analyzing Fluorescence Recovery After
%     Photobleaching (FRAP) data. Copyright (C) 2011  Lewis J. Kraft
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function Figure_GUI_FRAP_FRET(data, basicinput)
% This is a figure GUI for the diffusion model
%--------------------------------------------------------------------------
%% Format GUI items

GUIfigureh = figure('Position',[200 200 1000 400],...
    'MenuBar','none','NumberTitle','off',...
    'Name','Data Analysis and Visualization','Visible','off');

%--------------------------------------------------------------------------
UserInputsh = uipanel('Title','User Inputs','Units',...
    'normalized','Position',[0.0169    0.0169    0.43    0.97],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14,'FontWeight','bold');

uicontrol('Parent',UserInputsh,'Style','text',...
    'String','Define the curve fitting parameters:','Units','normalized','Position',[0.025    0.9    .975    0.0741],...
    'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'));

dat =  {1,  0,  Inf,    'Adjustable';...
    3,  0,  Inf,    'Adjustable';...
    10, 0,  Inf,    'Adjustable';
    1,  0,  2,  'Adjustable';
    -0.001, -Inf,   0,  'Adjustable'};
columnname = {'Initial Guess', 'Lower Bound', 'Upper Bound', 'Fixed/Adj'};
rowname =   {'K', 're', 'D', 'MF', 'kdecay'};
columnformat = {'numeric', 'numeric', 'numeric', {'Fixed' 'Adjustable'}};
columneditable =  [true true true true];
t1 = uitable('Parent',UserInputsh,'Units','normalized','Position',...
    [.025,.6,.95,.31], 'Data', dat,'ColumnName', columnname,...
    'ColumnFormat', columnformat,'ColumnEditable', columneditable,...
    'RowName',rowname);
%--------------------------------------------------------------------------
uicontrol('Parent',UserInputsh,'Style','text',...
    'String','Define data inclusion parameters:','Units','normalized','Position',[0.025    0.45    .975    0.0741],...
    'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'));

dat =  {1,  length(data(1).r);
    1,  length(data(1).frap);
    length(data(1).frap)-50,   length(data(1).frap);
    length(data(1).frap)-.1*length(data(1).frap), length(data(1).frap)};
columnname = {'First data point','Last data point'};
rowname =   {'Profile Fit', 'FRAP Fit', 'Decay Fit', 'Corrected MF'};
columnformat = {'numeric', 'numeric', 'numeric', 'numeric'};
columneditable =  [true true true true true];
t2 = uitable('Parent',UserInputsh,'Units','normalized','Position',...
    [.025,.25,.95,.21], 'Data', dat,'ColumnName', columnname,...
    'ColumnFormat', columnformat,'ColumnEditable', columneditable,...
    'RowName',rowname);
%--------------------------------------------------------------------------
% Normh = uicontrol('Style','edit','String',...
%     '20','Units','normalized','Position',[.2875,.3,.15,.08]);
%--------------------------------------------------------------------------
uicontrol('Parent',UserInputsh,'Style','text',...
    'String','Fit the averaged data?','Units','normalized','Position',[0.15    0.1    .85    0.06],...
    'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'),'HorizontalAlignment','left');
Avgh = uicontrol('Parent',UserInputsh,'Style','popupmenu','String',{'No','Yes'},'Units','normalized','Position',[.6,.1,.2,.06]);
%--------------------------------------------------------------------------
FitOutputsh = uipanel('Title','Fit Outputs','Units',...
    'normalized','Position',[0.46    0.0169    0.43    0.97],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14,'FontWeight','bold');

columnname = {'K', 're', 'D', 'MF', 'MF Correct', 'SS', 'E','FRET Slope', 'FRET SS'};
rowname =   [{basicinput{:,1}},{'Avg.'}];
columnformat = [repmat({'numeric'},1,8),{'numeric'}];
columneditable =  [false false false false false false false false false];
t3 = uitable('Parent',FitOutputsh,'Units','normalized','Position',...
    [.025,.025,.95,.95], ...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', columneditable,...
    'RowName',rowname);
%--------------------------------------------------------------------------
listboxh = uicontrol(GUIfigureh,'Style','listbox','Units','normalized',...
    'String',basicinput(:,1),'Max',2,'Min',0,'Position',[.9, .155, .075, .8]);
plotbuttonh = uicontrol(GUIfigureh,'Style','pushbutton','Units','normalized',...
    'String','Run','Position',[.9 .085 .075 .06],...
    'Callback',{@plotbutton_callback,data},'FontSize',14);
savebuttonh = uicontrol(GUIfigureh,'Style','pushbutton','Units','normalized',...
    'String','Save','Position',[.9 .025 .075 .06],...
    'Callback',{@savebutton_callback,data},'FontSize',14);

%--------------------------------------------------------------------------

%% Define the Callback functions

    function data=plotbutton_callback(hObj,event,data)
        
        pbpplotsh=figure('Name','Post-bleach profile plots','NumberTitle','off');
        frapplotsh=figure('Name','FRAP plots','NumberTitle','off');
        fretplotsh=figure('Name','FRET plots','NumberTitle','off');
        normfretplotsh=figure('Name','FRET plots','NumberTitle','off');
        %         cla(resFRAPploth); cla(avgFRAPploth); cla(resavgFRAPploth); cla(pbpplotsh); cla(frapplotsh);
        val=get(listboxh,'Value');
        linecolors=lines(length(val));
        usrinputs=get(t1,'Data');
        usrinputs(6:9,1:2)=get(t2,'Data');
        usrinputs{10,1}=get(Avgh,'Value');
        assignin('base', 'usrinputs', usrinputs);
        %------------------------------------------------------------------
        % Correcting for unintentional photobleaching----------------------
        switch usrinputs{5,4}
            case 'Adjustable'
                [decayrate data]=PhotoDecay_FRAP_FRET(data,basicinput,usrinputs,val);
            case 'Fixed' %If you don't want to correct for photodecay enter 0.
                for index1=1:length(val)
                    t=data(val(index1)).time-data(val(index1)).time(1)';
                    data(val(index1)).correctfrap=data(val(index1)).normfrap./exp(usrinputs{5,1}*t(1,:));
                    data(val(index1)).correctfret=data(val(index1)).normfret./exp(usrinputs{5,1}*t(1,:));
                    decayrate=usrinputs{5,1};
                end
        end
        assignin('base', 'data', data);
        assignin('base', 'decayrate', decayrate);
        % End Correcting for unintentional photobleaching------------------
        %------------------------------------------------------------------
        
        % Fit FRAP/FRET Data
        [data avg]=FRAP_FRET_Model(data,basicinput,usrinputs,val);
        
        
        assignin('base', 'dataout', data);
        assignin('base', 'avg', avg);
        %------------------------------------------------------------------
        figure(fretplotsh);
        subplot(3,2,[1,3])
        for index1=1:length(val)
            line(data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),data(val(index1)).correctfret(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),'Line','none','Marker','o','Color',linecolors(index1,:))
            if usrinputs{10,1}==1
                line(0:.1:1.1,data(val(index1)).fretfit,'Color','k','LineWidth',2)
            end
        end
        ylabel({'Donor Intensity','(normalized)'})
        grid on
        subplot(3,2,[5])
        if usrinputs{10,1}==1
            for index1=1:length(val)
                line(data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),data(val(index1)).fretres,'Line','none','Marker','o','Color',linecolors(index1,:))
            end
        end
        ylabel({'Residuals'})
        xlabel({'Acceptor Intensity (normalized)'})
        grid on
        subplot(3,2,[2,4])
        line(avg.f(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),avg.f2(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),'Line','none','Marker','o','Color','r')
        line(0:.1:1.1,avg.fretfit,'Color','k','LineWidth',2)
        grid on
        subplot(3,2,[6])
        line(avg.f(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}),avg.fretres,'Line','none','Marker','o','Color',linecolors(index1,:))
        xlabel({'Acceptor Intenisty (normalized)'})
        grid on
        %------------------------------------------------------------------
        figure(normfretplotsh);
        subplot(1,2,1)
        for index1=1:length(val)
            line(data(val(index1)).time-data(val(index1)).time(basicinput{1,6}),data(val(index1)).correctfret,'Line','none','Marker','o','Color',linecolors(index1,:))
        end
        ylabel({'Donor Fluorescence Intensity','(normalized)'})
        xlabel({'Time (s)'})
        grid on
        subplot(1,2,2)
        line(avg.time-avg.time(basicinput{1,6}),avg.f2,'Line','none','Marker','o','Color','r')
        xlabel({'Time (s)'})
        grid on
        %------------------------------------------------------------------
        
        figure(frapplotsh);
        subplot(3,2,[1,3])
        for index1=1:length(val)
            line(data(val(index1)).time-data(val(index1)).time(basicinput{1,6}),data(val(index1)).correctfrap,'Line','none','Marker','o','Color',linecolors(index1,:))
            if usrinputs{10,1}==1
                line(data(val(index1)).t,data(val(index1)).frapfit,'Color','k','LineWidth',2)
            end
        end
        ylabel({'Fluorescence Intensity','(normalized)'})
        grid on
        subplot(3,2,[5])
        if usrinputs{10,1}==1
            for index1=1:length(val)
                line(data(val(index1)).t,data(val(index1)).frapres,'Line','none','Marker','o','Color',linecolors(index1,:))
            end
        end
        ylabel({'Residuals'})
        xlabel({'Time (s)'})
        grid on
        subplot(3,2,[2,4])
        line(avg.time-avg.time(basicinput{1,6}),avg.f,'Line','none','Marker','o','Color','r')
        line(avg.t,avg.frapfit,'Color','k','LineWidth',2)
        grid on
        subplot(3,2,[6])
        line(avg.t,avg.frapres,'Line','none','Marker','o','Color',linecolors(index1,:))
        xlabel({'Time (s)'})
        grid on
        %------------------------------------------------------------------
        
        figure(pbpplotsh)
        subplot(3,2,[1,3])
        for index1=1:length(val)
            line(data(val(index1)).r,data(val(index1)).pbp,'Line','none','Marker','o','Color',linecolors(index1,:))
            if usrinputs{10,1}==1
                line(data(val(index1)).rfit,data(val(index1)).pbpfit,'Color','k','LineWidth',2)
            end
        end
        ylabel({'Fluorescence Intensity','(normalized)'})
        grid on
        subplot(3,2,[5])
        if usrinputs{10,1}==1
            for index1=1:length(val)
                line(data(val(index1)).rfit,data(val(index1)).pbpres,'Line','none','Marker','o','Color',linecolors(index1,:))
            end
        end
        ylabel({'Residuals'})
        xlabel({'Radial distance (\mum)'})
        grid on
        subplot(3,2,[2,4])
        plot(avg.r,avg.pbp,'.k',avg.rfit,avg.pbpfit,'r','LineWidth',2)
        grid on
        subplot(3,2,[6])
        plot(avg.rfit,avg.pbpres,'.k')
        xlabel('Radial distance (\mum)')
        grid on
        %------------------------------------------------------------------
        % upload the fit parameters into table2
        if usrinputs{10,1}==1
            for index1=1:length(val)
                temp{val(index1),1}=data(val(index1)).k;
                temp{val(index1),2}=data(val(index1)).re;
                temp{val(index1),3}=data(val(index1)).D;
                temp{val(index1),4}=data(val(index1)).MF;
                temp{val(index1),5}=data(val(index1)).correctMF;
                temp{val(index1),6}=data(val(index1)).SS;
                temp{val(index1),7}=data(val(index1)).E;
                temp{val(index1),8}=data(val(index1)).fretslope;
                temp{val(index1),9}=data(val(index1)).fretSS;
            end
        end
        temp{length([basicinput(:,1)])+1,1}=avg.k;
        temp{length([basicinput(:,1)])+1,2}=avg.re;
        temp{length([basicinput(:,1)])+1,3}=avg.D;
        temp{length([basicinput(:,1)])+1,4}=avg.MF;
        temp{length([basicinput(:,1)])+1,5}=avg.correctMF;
        temp{length([basicinput(:,1)])+1,6}=avg.SS;
        temp{length([basicinput(:,1)])+1,7}=avg.E;
        temp{length([basicinput(:,1)])+1,8}=avg.fretslope;
        temp{length([basicinput(:,1)])+1,9}=avg.fretSS;
        
        set(t3,'Data',temp);
        assignin('base','temp',temp);
        clearvars temp
        %------------------------------------------------------------------
        
    end

    function data=savebutton_callback(hObj,event,data)
        answer = inputdlg({'Enter a file description:'},'Save File Descriptor',1,{'Venus_Cytoplasm'});
        temp=evalin('base','temp');
        rowname=get(t3,'RowName');
        assignin('base','rowname',rowname);
        usrinputs=evalin('base', 'usrinputs');
        FileLocation=evalin('base','FileLocation');
        savedata=[rowname, temp]';
        header={'FileNames','k','re','D','MF','MF Corrected','SS','E','FRET Slope','FRET SS'};
        fid = fopen(fullfile(FileLocation,[answer{:},'_FRAP_FRET_Fit_Parameters.txt']),'w');
        fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r\n', header{:});
        fprintf(fid, '%s\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\r\n', savedata{:});
        fclose(fid);
        
        dataout=evalin('base','dataout');
        avg=evalin('base','avg');
        
        m=length(dataout)+1;
        count=1;
        %if exist dataout.t  % this needs to be inserted to prevent an
        %error in the fitting of average data.
        for j=1:length(dataout)
            a{count}=sprintf([basicinput{j,1},'\t Time\t',repmat('%f\t',1,length(dataout(j).time)-1),'%f\r\n'],dataout(j).time-dataout(j).time(basicinput{1,6}));
            a{count+1}=sprintf([basicinput{j,1},'\t Raw FRAP\t',repmat('%f\t',1,length(dataout(j).frap)-1),'%f\r\n'],dataout(j).frap);
            a{count+2}=sprintf([basicinput{j,1},'\t Cell\t',repmat('%f\t',1,length(dataout(j).cell)-1),'%f\r\n'],dataout(j).cell);
            a{count+3}=sprintf([basicinput{j,1},'\t Normalized FRAP\t',repmat('%f\t',1,length(dataout(j).normfrap)-1),'%f\r\n'],dataout(j).normfrap);
            a{count+4}=sprintf([basicinput{j,1},'\t Corrected FRAP\t',repmat('%f\t',1,length(dataout(j).correctfrap)-1),'%f\r\n'],dataout(j).correctfrap);
            a{count+5}=sprintf([basicinput{j,1},'\t Fit Time\t',repmat('%f\t',1,length(dataout(j).t)-1),'%f\r\n'],dataout(j).t);
            a{count+6}=sprintf([basicinput{j,1},'\t FRAP Fit\t',repmat('%f\t',1,length(dataout(j).frapfit)-1),'%f\r\n'],dataout(j).frapfit);
            a{count+7}=sprintf([basicinput{j,1},'\t Fit Residuals\t',repmat('%f\t',1,length(dataout(j).frapres)-1),'%f\r\n'],dataout(j).frapres);
            a{count+8}=sprintf([basicinput{j,1},'\t Acceptor Intensity\t',repmat('%f\t',1,length(dataout(j).correctfrap(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}))-1),'%f\r\n'],dataout(j).correctfrap(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}));
            a{count+9}=sprintf([basicinput{j,1},'\t Donor Intensity\t',repmat('%f\t',1,length(dataout(j).correctfret(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}))-1),'%f\r\n'],dataout(j).correctfret(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2}));
            a{count+10}=sprintf([basicinput{j,1},'\t FRET Fit Time\t',repmat('%f\t',1,length([0:.1:1.1])-1),'%f\r\n'],[0:.1:1.1]);
            a{count+11}=sprintf([basicinput{j,1},'\t FRET Fit\t',repmat('%f\t',1,length(dataout(j).fretfit)-1),'%f\r\n'],dataout(j).fretfit);
            a{count+12}=sprintf([basicinput{j,1},'\t FRET Fit Residuals\t',repmat('%f\t',1,length(dataout(j).fretres)-1),'%f\r\n'],dataout(j).fretres);
            count=count+13;
        end
        j=1;
        a{count}=sprintf(['Average\t Time\t',repmat('%f\t',1,length(avg(j).time)-1),'%f\r\n'],avg(j).time-avg(j).time(basicinput{1,6}));
        a{count+1}=sprintf(['Average\t Corrected FRAP\t',repmat('%f\t',1,length(avg(j).f)-1),'%f\r\n'],avg(j).f);
        a{count+2}=sprintf(['Average\t Fit Time\t',repmat('%f\t',1,length(avg(j).t)-1),'%f\r\n'],avg(j).t);
        a{count+3}=sprintf(['Average\t FRAP Fit\t',repmat('%f\t',1,length(avg(j).frapfit)-1),'%f\r\n'],avg(j).frapfit);
        a{count+4}=sprintf(['Average\t Fit Residuals\t',repmat('%f\t',1,length(avg(j).frapres)-1),'%f\r\n'],avg(j).frapres);
        a{count+5}=sprintf(['Average\t Acceptor Intensity\t',repmat('%f\t',1,length(avg(j).f)-1),'%f\r\n'],avg(j).f);
        a{count+6}=sprintf(['Average\t Donor Intensity\t',repmat('%f\t',1,length(avg(j).f2)-1),'%f\r\n'],avg(j).f2);
        a{count+7}=sprintf(['Average\t FRET Fit Time\t',repmat('%f\t',1,length([0:.1:1.1])-1),'%f\r\n'],[0:.1:1.1]);
        a{count+8}=sprintf(['Average\t FRET Fit\t',repmat('%f\t',1,length(avg(j).fretfit)-1),'%f\r\n'],avg(j).fretfit);
        a{count+9}=sprintf(['Average\t FRET Fit Residuals\t',repmat('%f\t',1,length(avg(j).fretres)-1),'%f\r\n'],avg(j).fretres);
            
        str=[a{:}];
        fid = fopen(fullfile(FileLocation,[answer{:},'_FRAP_FRET_datasets.txt']),'w');
        fprintf(fid,str);
        fclose(fid);
        clearvars a str
        
        m=length(dataout)+1;
        count=1;
        for j=1:length(dataout)
            a{count}=sprintf([basicinput{j,1},'\t Distance\t',repmat('%f\t',1,length(dataout(j).r)-1),'%f\r\n'],dataout(j).r);
            a{count+1}=sprintf([basicinput{j,1},'\t Post-bleach Profile\t',repmat('%f\t',1,length(dataout(j).pbp)-1),'%f\r\n'],dataout(j).pbp);
            a{count+2}=sprintf([basicinput{j,1},'\t Profile Fit\t',repmat('%f\t',1,length(dataout(j).pbpfit)-1),'%f\r\n'],dataout(j).pbpfit);
            a{count+3}=sprintf([basicinput{j,1},'\t Fit Residuals\t',repmat('%f\t',1,length(dataout(j).pbpres)-1),'%f\r\n'],dataout(j).pbpres);
            count=count+4;
        end
        j=1;
        a{count}=sprintf(['Average\t Distance\t',repmat('%f\t',1,length(avg(j).r)-1),'%f\r\n'],avg(j).r);
        a{count+1}=sprintf(['Average\t Post-bleach Profile\t',repmat('%f\t',1,length(avg(j).pbp)-1),'%f\r\n'],avg(j).pbp);
        a{count+2}=sprintf(['Average\t Profile Fit\t',repmat('%f\t',1,length(avg(j).pbpfit)-1),'%f\r\n'],avg(j).pbpfit);
        a{count+3}=sprintf(['Average\t Fit Residuals\t',repmat('%f\t',1,length(avg(j).pbpres)-1),'%f\r\n'],avg(j).pbpres);
        
        str=[a{:}];
        fid = fopen(fullfile(FileLocation,[answer{:},'_FRAP_FRET_Postbleach_profiles.txt']),'w');
        fprintf(fid,str);
        fclose(fid);
        clearvars a str
    end

% Initialize the GUI.
% Change units to normalized so components resize
% automatically.
%--------------------------------------------------------------------------
% Move the GUI to the center of the screen.
movegui(GUIfigureh,'center')
% Make the GUI visible.
set(GUIfigureh,'Visible','on');

%--------------------------------------------------------------------------

end