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

function Figure_GUI_Reaction2(data, basicinput)
% This is a figure GUI for the Reaction 2 model
% Inputs:
% data - this is the data output from loadData_Reaction.m it contains the
% processed FRAP curves.
% basicinput - basic user input from main_GUI.

% Outputs:
%

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

dat =  {1,  0,  Inf,    'Adjustable';... % These are default values for inital, lower, and upper bounds.
    0.5,  0,  Inf,    'Adjustable';
    0.05, 0, Inf, 'Adjustable';
    0.5, 0, Inf, 'Adjustable';
    0.0005, 0, Inf, 'Adjustable';
    0.001, 0, Inf, 'Adjustable'};

columnname = {'Initial Guess', 'Lower Bound', 'Upper Bound', 'Fixed/Adj'};
rowname =   {'a', 'b', 'c', 'd', 'f', 'decayrate'};
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

dat =  { 1 , length(data(1).frap) ; %These are the default values for data exclusion
    length(data(1).frap)-50,   length(data(1).frap)};

columnname = {'First data point','Last data point'};
rowname =   {'FRAP Fit', 'Decay Fit'};
columnformat = {'numeric', 'numeric'};
columneditable =  [true true];
t2 = uitable('Parent',UserInputsh,'Units','normalized','Position',...
    [.025,.25,.95,.21], 'Data', dat,'ColumnName', columnname,...
    'ColumnFormat', columnformat,'ColumnEditable', columneditable,...
    'RowName',rowname);
%--------------------------------------------------------------------------
% Normh = uicontrol('Style','edit','String',...
%     '20','Units','normalized','Position',[.2875,.3,.15,.08]);
%--------------------------------------------------------------------------
uicontrol('Parent',UserInputsh,'Style','text',... % This prompts the user if they would like to fit the averaged data or the individual data
    'String','Fit the averaged data?','Units','normalized','Position',[0.15    0.1    .85    0.06],...
    'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'),'HorizontalAlignment','left');
Avgh = uicontrol('Parent',UserInputsh,'Style','popupmenu','String',{'No','Yes'},'Units','normalized','Position',[.6,.1,.2,.06]);
%--------------------------------------------------------------------------
FitOutputsh = uipanel('Title','Fit Outputs','Units',...
    'normalized','Position',[0.46    0.0169    0.43    0.97],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14,'FontWeight','bold');

columnname = {'a', 'b', 'c', 'd', 'f', 'SS'}; % This is the table where the optimized parameters will output for visual inspection
rowname =   [{basicinput{:,1}},{'Avg.'}];
columnformat = {'numeric','numeric','numeric','numeric','numeric','numeric'};
columneditable =  [false false false false false false];
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
        
        frapplotsh=figure('Name','FRAP plots','NumberTitle','off'); % Create a new figure window for visualizing the processed data and the optimized fits.
        val=get(listboxh,'Value'); % This fetches the particular datasets that the user selected to plot
        linecolors=lines(length(val)); % Each frapcurve will be plotted with a different color.
        usrinputs=get(t1,'Data'); % Fetch the initial, lower, and upper bounds.
        usrinputs(7:8,1:2)=get(t2,'Data'); % Fetch data exclusion parameters.
        usrinputs{9,1}=get(Avgh,'Value'); % Fetches the choice about fitting the average data
        assignin('base', 'usrinputs', usrinputs);
        
        %------------------------------------------------------------------
        %% Correcting for unintentional photobleaching----------------------
        % there should be a conditional to not do this if we normalized by
        % whole cell
        if basicinput{1,4}==2
            for index1=1:length(val)
                data(val(index1)).correctfrap=data(val(index1)).normfrap;
                decayrate=0;
            end
        else
        switch usrinputs{6,4}
            case 'Adjustable'
                [decayrate data]=PhotoDecay_Reaction2(data,basicinput,usrinputs,val); % Correct the photodecay in the FRAP datasets.
            case 'Fixed' %If you don't want to correct for photodecay select 'Fixed' in the initial, lower, and upper bound box, and enter 0 as the initial value.
                for index1=1:length(val)
                    t=data(val(index1)).time-data(val(index1)).time(1)';
                    data(val(index1)).correctfrap=data(val(index1)).normfrap./exp(-usrinputs{6,1}*t(1,:));
                    decayrate=usrinputs{6,1};
                end
        end
        end
        
        
        assignin('base', 'data', data);
        assignin('base', 'decayrate', decayrate);
        % End Correcting for unintentional photobleaching------------------
        
        %------------------------------------------------------------------
        
        %% Fit with one component reaction model
        [data, avg]=ReactionModel2(data,basicinput,usrinputs,val);
        assignin('base', 'dataout', data);
        assignin('base', 'avg', avg);
        %------------------------------------------------------------------
        
        %% Visualize the results of the fitting
        figure(frapplotsh);
        % Plot all of the individual FRAP datasets and optimized fits
        subplot(3,2,[1,3])
        for index1=1:length(val)
            line(data(val(index1)).time-data(val(index1)).time(basicinput{1,6}),data(val(index1)).correctfrap,'Line','none','Marker','o','Color',linecolors(index1,:))
            if usrinputs{9,1}==1
                line(data(val(index1)).t,data(val(index1)).frapfit,'Color','k','LineWidth',2)
            end
        end
        ylabel({'Fluorescence Intensity','(normalized)'})
        grid on
        subplot(3,2,[5])
        if usrinputs{9,1}==1
            for index1=1:length(val)
                line(data(val(index1)).t,data(val(index1)).frapres,'Line','none','Marker','o','Color',linecolors(index1,:))
            end
        end
        ylabel({'Residuals'})
        xlabel({'Time (s)'})
        grid on
        % Plot the average of all FRAP datasets and the averaged optimized
        % fit
        subplot(3,2,[2,4])
        line(avg.time-avg.time(basicinput{1,6}),avg.f,'Line','none','Marker','o','Color','r')
        line(avg.t,avg.frapfit,'Color','k','LineWidth',2)
        grid on
        subplot(3,2,[6])
        line(avg.t,avg.frapres,'Line','none','Marker','o','Color',linecolors(index1,:))
        xlabel({'Time (s)'})
        grid on
        %------------------------------------------------------------------
        
        %% upload the fit parameters into t3 (table 3)
        if usrinputs{9,1}==1
            for index1=1:length(val)
                temp{val(index1),1}=data(val(index1)).a;
                temp{val(index1),2}=data(val(index1)).b;
                temp{val(index1),3}=data(val(index1)).c;
                temp{val(index1),4}=data(val(index1)).d;
                temp{val(index1),5}=data(val(index1)).fparam;
                temp{val(index1),6}=data(val(index1)).SS;
            end
        end
        temp{length([basicinput(:,1)])+1,1}=avg.a;
        temp{length([basicinput(:,1)])+1,2}=avg.b;
        temp{length([basicinput(:,1)])+1,3}=avg.c;
        temp{length([basicinput(:,1)])+1,4}=avg.d;
        temp{length([basicinput(:,1)])+1,5}=avg.fparam;
        temp{length([basicinput(:,1)])+1,6}=avg.SS;
        
        set(t3,'Data',temp);
        assignin('base','temp',temp);
        clearvars temp
        %------------------------------------------------------------------
        
    end


%% Saving the optimized parameters, FRAP curves, and FRAP fits to a tab delimited text file
    function data=savebutton_callback(hObj,event,data)
        answer = inputdlg({'Enter a file description:'},'Save File Descriptor',1,{'Venus_Cytoplasm'});
        temp=evalin('base','temp');
        rowname=get(t3,'RowName'); % Fetch the optimized parameters from table 3
        assignin('base','rowname',rowname);
        FileLocation=evalin('base','FileLocation');
        savedata=[rowname, temp]';
        header={'FileNames','a','b','c','d','f','SS'};
        fid = fopen(fullfile(FileLocation,[answer{:},'_Reaction2_Fit_Parameters.txt']),'w');
        fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\r\n', header{:});
        fprintf(fid, '%s\t %g\t %g\t %g\t %g\t %g\t %g\r\n', savedata{:});
        fclose(fid);
        
        dataout=evalin('base','dataout');
        avg=evalin('base','avg');
        
        m=length(dataout)+1;
        count=1;
        for j=1:length(dataout)
            a{count}=sprintf([basicinput{j,1},'\t Time\t',repmat('%f\t',1,length(dataout(j).time)-1),'%f\r\n'],dataout(j).time-dataout(j).time(basicinput{1,6}));
            a{count+1}=sprintf([basicinput{j,1},'\t Raw FRAP\t',repmat('%f\t',1,length(dataout(j).frap)-1),'%f\r\n'],dataout(j).frap);
            a{count+2}=sprintf([basicinput{j,1},'\t Cell\t',repmat('%f\t',1,length(dataout(j).cell)-1),'%f\r\n'],dataout(j).cell);
            a{count+3}=sprintf([basicinput{j,1},'\t Normalized FRAP\t',repmat('%f\t',1,length(dataout(j).normfrap)-1),'%f\r\n'],dataout(j).normfrap);
            a{count+4}=sprintf([basicinput{j,1},'\t Corrected FRAP\t',repmat('%f\t',1,length(dataout(j).correctfrap)-1),'%f\r\n'],dataout(j).correctfrap);
            a{count+5}=sprintf([basicinput{j,1},'\t Fit Time\t',repmat('%f\t',1,length(dataout(j).t)-1),'%f\r\n'],dataout(j).t);
            a{count+6}=sprintf([basicinput{j,1},'\t FRAP Fit\t',repmat('%f\t',1,length(dataout(j).frapfit)-1),'%f\r\n'],dataout(j).frapfit);
            a{count+7}=sprintf([basicinput{j,1},'\t Fit Residuals\t',repmat('%f\t',1,length(dataout(j).frapres)-1),'%f\r\n'],dataout(j).frapres);
            count=count+8;
        end
        j=1;
        a{count}=sprintf(['Average\t Time\t',repmat('%f\t',1,length(avg(j).time)-1),'%f\r\n'],avg(j).time-avg(j).time(basicinput{1,6}));
        a{count+1}=sprintf(['Average\t Corrected FRAP\t',repmat('%f\t',1,length(avg(j).f)-1),'%f\r\n'],avg(j).f);
        a{count+2}=sprintf(['Average\t Fit Time\t',repmat('%f\t',1,length(avg(j).t)-1),'%f\r\n'],avg(j).t);
        a{count+3}=sprintf(['Average\t FRAP Fit\t',repmat('%f\t',1,length(avg(j).frapfit)-1),'%f\r\n'],avg(j).frapfit);
        a{count+4}=sprintf(['Average\t Fit Residuals\t',repmat('%f\t',1,length(avg(j).frapres)-1),'%f\r\n'],avg(j).frapres);
        
        str=[a{:}];
        fid = fopen(fullfile(FileLocation,[answer{:},'_Reaction2_FRAP_datasets.txt']),'w');
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