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

function Figure_GUI_Diffusion(data, basicinput)
% This is a figure GUI for the diffusion model
% Inputs:
% data - this is the output from loadData_Diffusion.m
% basicinput - basic user input from the Main_GUI

% Outputs:
%

%--------------------------------------------------------------------------
%% Format GUI items
screen_size = get(0, 'ScreenSize');
GUIfigureh = figure('Position',[0,0,screen_size(3),screen_size(4)],...
    'MenuBar','none','NumberTitle','off',...
    'Name','Data Analysis and Visualization','Visible','off');

%--------------------------------------------------------------------------
UserInputsh = uipanel('Title','User Inputs','Units',...
    'normalized','Position',[0.0169    0.0169    0.32    0.97],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14,'FontWeight','bold');

uicontrol('Parent',UserInputsh,'Style','text',...
    'String','Define the curve fitting parameters:','Units','normalized','Position',[0.025    0.9    .975    0.0741],...
    'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'));

dat =  {1,  0,  Inf,    'Adjustable';... % These are the default intial, lower, and upper bounds on the fitting parameters.
    3,  0,  Inf,    'Adjustable';...
    10, 0,  Inf,    'Adjustable';
    1,  0,  2,  'Adjustable';
    0.001, 0,   Inf,  'Adjustable'};
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

for j=1:length(data);
    rSize(j)=length(data(j).r);
end
dat =  {1,  min(rSize); % These are the default data exclusion parameters.
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
    'normalized','Position',[0.35    0.0169    0.52    0.97],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14,'FontWeight','bold');

columnname = {'K', 're', 'D', 'MF', 'MF Correct', 'SS'}; % The optimized parameters get uploaded into this table
rowname =   [{basicinput{:,1}},{'Avg.'}];
columnformat = [repmat({'numeric'},1,5),{'numeric'}];
columneditable =  [false false false false false false];
t3 = uitable('Parent',FitOutputsh,'Units','normalized','Position',...
    [.025,.025,.95,.95], ...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', columneditable,...
    'RowName',rowname);
%--------------------------------------------------------------------------
listboxh = uicontrol(GUIfigureh,'Style','listbox','Units','normalized',...
    'String',basicinput(:,1),'Max',2,'Min',0,'Position',[.88, .155, .1, .8]);
plotbuttonh = uicontrol(GUIfigureh,'Style','pushbutton','Units','normalized',...
    'String','Run','Position',[.9 .085 .075 .06],...
    'Callback',{@plotbutton_callback,data},'FontSize',14);
savebuttonh = uicontrol(GUIfigureh,'Style','pushbutton','Units','normalized',...
    'String','Save','Position',[.9 .025 .075 .06],...
    'Callback',{@savebutton_callback,data},'FontSize',14);

%--------------------------------------------------------------------------

%% Define the Callback functions

    function data=plotbutton_callback(hObj,event,data)
        %%
        try
            pbpplotsh=figure('Name','Post-bleach profile plots','NumberTitle','off'); % Create figure window for plotting the intial conditions.
            frapplotsh=figure('Name','FRAP plots','NumberTitle','off'); % Create figure window for plotting the FRAP curves.
            val=get(listboxh,'Value'); % Fetch the FRAP datasets selected by the user for analysis.
            linecolors=lines(length(val)); % Plot each dataset with its own unique color.
            usrinputs=get(t1,'Data'); % Fetch the initial, lower, and upper bounds from table 1.
            usrinputs(6:9,1:2)=get(t2,'Data'); % Fetch the data exclusion parameters from table 2.
            usrinputs{10,1}=get(Avgh,'Value'); % The user decided if they want to fit all datasets individually or average them before fitting.
            assignin('base', 'usrinputs', usrinputs);
            %------------------------------------------------------------------
            
            
            %% Correcting for unintentional photobleaching----------------------
            if basicinput{1,4}==2
                for index1=1:length(val)
                    data(val(index1)).correctfrap=data(val(index1)).normfrap;
                    decayrate=0;
                end
            else
                switch usrinputs{5,4}
                    case 'Adjustable' % If the user wants to fit the end of the FRAP curve to find the rate of photodecay.
                        [decayrate data]=PhotoDecay(data,basicinput,usrinputs,val); % Correct the photodecay.
                    case 'Fixed' %If you don't want to correct for photodecay fix the photodecay rate and enter 0.
                        for index1=1:length(val)
                            t=data(val(index1)).time-data(val(index1)).time(1)';
                            data(val(index1)).correctfrap=data(val(index1)).normfrap./exp(-usrinputs{5,1}*t(1,:));
                            decayrate=usrinputs{5,1};
                        end
                end
            end
            
            
            assignin('base', 'data', data);
            assignin('base', 'decayrate', decayrate);
            % End Correcting for unintentional photobleaching------------------
            %------------------------------------------------------------------
            
            %% Fit with diffusion model from Kang et al, Biophys J, 2008.
            [data avg]=DiffusionModel_2(data,basicinput,usrinputs,val);
            
            
            assignin('base', 'dataout', data);
            assignin('base', 'avg', avg);
            %------------------------------------------------------------------
            
            %% Plot the results for visualization
            figure(frapplotsh); % Plot the FRAP curves
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
            
            figure(pbpplotsh) % Plot the initial conditions
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
            
            
            %% upload the optimized fit parameters into table 3
            if usrinputs{10,1}==1
                for index1=1:length(val)
                    temp{val(index1),1}=data(val(index1)).k;
                    temp{val(index1),2}=data(val(index1)).re;
                    temp{val(index1),3}=data(val(index1)).D;
                    temp{val(index1),4}=data(val(index1)).MF;
                    temp{val(index1),5}=data(val(index1)).correctMF;
                    temp{val(index1),6}=data(val(index1)).SS;
                end
            end
            temp{length([basicinput(:,1)])+1,1}=avg.k;
            temp{length([basicinput(:,1)])+1,2}=avg.re;
            temp{length([basicinput(:,1)])+1,3}=avg.D;
            temp{length([basicinput(:,1)])+1,4}=avg.MF;
            temp{length([basicinput(:,1)])+1,5}=avg.correctMF;
            temp{length([basicinput(:,1)])+1,6}=avg.SS;
            
            set(t3,'Data',temp);
            assignin('base','temp',temp);
            clearvars temp
            %------------------------------------------------------------------
        catch errObj
            errordlg(getReport(errObj,'extended','hyperlinks','off'),'Error 1');
            close(frapplotsh);
            close(pbpplotsh);
            return
        end
        
    end

%% Save the FRAP curves, optimized fit parameters, and initial conditions to a tab delimited text file.

    function data=savebutton_callback(hObj,event,data)
        try
            answer = inputdlg({'Enter a file description:'},'Save File Descriptor',1,{'Venus_Cytoplasm'});
            temp=evalin('base','temp');
            rowname=get(t3,'RowName');
            assignin('base','rowname',rowname);
            FileLocation=evalin('base','FileLocation');
            savedata=[rowname, temp]';
            header={'FileNames','k','re','D','MF','MF Corrected','SS'};
            fid = fopen(fullfile(FileLocation,[answer{:},'_Diffusion_Fit_Parameters.txt']),'w');
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
            fid = fopen(fullfile(FileLocation,[answer{:},'_Diffusion_FRAP_datasets.txt']),'w');
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
            fid = fopen(fullfile(FileLocation,[answer{:},'_Diffusion_Postbleach_profiles.txt']),'w');
            fprintf(fid,str);
            fclose(fid);
            clearvars a str
            
        catch errObj
            errordlg(getReport(errObj,'extended','hyperlinks','off'),'Error 1');
            return
        end
        end
        
        %% Initialize the GUI.
        % Change units to normalized so components resize
        % automatically.
        %--------------------------------------------------------------------------
        % Move the GUI to the center of the screen.
        movegui(GUIfigureh,'center')
        % Make the GUI visible.
        set(GUIfigureh,'Visible','on');
        
        %--------------------------------------------------------------------------
        
    end