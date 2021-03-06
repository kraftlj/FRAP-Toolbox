%%% FRAP Toolbox
%%
%     FRAP Toolbox is designed to be a modular software program
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
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Running the script launches the Main_GUI - the beginning of FRAP-Toolbox.
%
function Main_GUI
% The user will be prompted for basic information about the FRAP
% experiments and the FRAP model to be used in downstream analyses.
%
%% Format the GUI
screen_size = get(0, 'ScreenSize');
GUIfigureh = figure('Position',[0,0,screen_size(3),screen_size(4)],...
    'MenuBar','none','NumberTitle','off',...
    'Name','FRAP Toolbox','Visible','off',...
    'Color',[.85,.85,.85]); % This creates a GUI window that is 800 x 400 pixels in size
%% Title and subtitle

uicontrol('Parent',GUIfigureh,'Style','text',...
    'String','FRAP Toolbox','Units','normalized','Position',[0.3912    0.9012    0.2189    0.0741],...
    'FontSize',24,'BackgroundColor',get(GUIfigureh,'color'),'ForegroundColor',[0,0,0],'FontWeight','bold');

%% Choose Directories Panel
ChooseDirectoriesh = uipanel('Title','Directory Information','Units',...
    'normalized','Position',[0.0169    0.6887    0.6587    0.2125],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14);

uicontrol('Parent',ChooseDirectoriesh,'Style','text',...
    'String','Choose directory containing FRAP data:','Units','normalized',...
    'Position',[0.025    .7    .6    .25],'HorizontalAlignment','left',...
    'BackgroundColor',get(GUIfigureh,'color'),'FontSize',14);

Filelocationh = uicontrol('Parent',ChooseDirectoriesh,'Style','edit',...
    'String','F:\Documents\Kenworthy','Units','normalized','Position',...
    [0.025    .25    .8    .3],'Callback',{@Filelocation_Callback},'FontSize',14);

Filelocationbuttonh = uicontrol('Parent',ChooseDirectoriesh,'Style',...
    'pushbutton','String','Browse','Units','normalized',...
    'Position',[.825,.25,.15,.3],'Callback',...
    {@Filelocationbutton_Callback},'FontSize',14);

    function Filelocationbutton_Callback(source,eventdata)
        set(FRAPfileslistboxh,'Value',[]);
        FRAP.FileLocation=uigetdir;
        set(Filelocationh,'String',FRAP.FileLocation);
        list=dir(fullfile(FRAP.FileLocation,'*'));
        set(FRAPfileslistboxh,'String',{list.name});
    end


%% User Inputs Panel

Userinputsh = uipanel('Title','User Inputs','Units','normalized',...
    'Position',[0.0169 .1 0.6587 .575],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14);

uicontrol('Parent',Userinputsh,'Style','text','String',...
    {'Select a FRAP model from the drop down list............................'},'Units','normalized','Position',...
    [.025,.91,.8,.08],'HorizontalAlignment','left','BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14);
FRAPmodelh = uicontrol('Parent',Userinputsh,'Style','popupmenu','String',...
    {'Diffusion','Reaction 1','Reaction 2','FRAP FRET'},'Units','normalized','Position',[.825,.91,.15,.08],'FontSize',14,'HorizontalAlignment','left');

uicontrol('Parent',Userinputsh,'Style','text',...
    'String',{'Select an ROI from the drop down list.........................................'},'Units',...
    'normalized','Position',[.025,.75,.8,.08],'HorizontalAlignment',...
    'left','BackgroundColor',get(GUIfigureh,'color'),'FontSize',14);
ROIh = uicontrol('Parent',Userinputsh,'Style','popupmenu','String',...
    {'Circle','User Defined'},'Units','normalized','Position',[.825,.75,.15,.08],'FontSize',14,'HorizontalAlignment','left');

uicontrol('Parent',Userinputsh,'Style','text','String',...
    {'Enter the frame number of the post-bleach image......................',...
    '(Valid entries are 2,3,4,...'},'Units','normalized',...
    'Position',[.025,.5,.8,.18],'HorizontalAlignment','left',...
    'BackgroundColor',get(GUIfigureh,'color'),'FontSize',14);
PBInumh = uicontrol('Parent',Userinputsh,'Style','edit','String','21',...
    'Units','normalized','Position',[.825,.59,.15,.08],'FontSize',14,'HorizontalAlignment','left');

uicontrol('Parent',Userinputsh,'Style','text','String'...
    ,{'Enter the mean background fluorescence intensity....................',...
    '(Valid entries are 0,1,2,...)'},'Units','normalized','Position',...
    [.025,.31,.8,.18],'HorizontalAlignment','left','BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14);
BGh = uicontrol('Parent',Userinputsh,'Style','edit','String',...
    '0','Units','normalized','Position',[.825,.4,.15,.08],'FontSize',14,'HorizontalAlignment','left');

uicontrol('Parent',Userinputsh,'Style','text','String',...
    {'Normalize by the fluorescence intensity of the whole cell?........'}...
    ,'Units','normalized','Position',[.025,.24,.8,.08],'HorizontalAlignment','left',...
    'BackgroundColor',get(GUIfigureh,'color'),'FontSize',14);
normalizationh = uicontrol('Parent',Userinputsh,'Style','popupmenu','String',...
    {'No','Yes'},'Units','normalized','Position',[.825,.24,.15,.08],'FontSize',14,'HorizontalAlignment','left');

uicontrol('Parent',Userinputsh,'Style','text','String',...
    {'How many pre-bleach images should be used for......................','normalization?'}...
    ,'Units','normalized','Position',[.025,.01,.8,.16],'HorizontalAlignment','left',...
    'BackgroundColor',get(GUIfigureh,'color'),'FontSize',14);
PBNh = uicontrol('Parent',Userinputsh,'Style','edit','String',...
    '10','Units','normalized','Position',[.825,.08,.15,.08],'FontSize',14,'HorizontalAlignment','left');

%% FRAP files in directory panel

DisplayFRAPfilesh = uipanel('Title','Select FRAP file(s)','Units',...
    'normalized','Position',[.6925 .1 .2906 .8012],'BackgroundColor',...
    get(GUIfigureh,'color'),'FontSize',14);

FRAPfileslistboxh = uicontrol('Parent',DisplayFRAPfilesh,'Style','listbox',...
    'String',{'Waiting for directory choice'},'Max',2,'Min',0,'Units',...
    'normalized','Position',[.025,.025,.95,.95],'FontSize',14);
    function Filelocation_Callback(source,eventdata)
        FRAP.FileLocation=get(Filelocationh,'String');
        list=dir(fullfile(FRAP.FileLocation,'*'));
        set(FRAPfileslistboxh,'String',{list.name});
    end


%% Analysis Buttons

Nextbuttonh = uicontrol('Parent',GUIfigureh,'Style','pushbutton',...
    'String','Next','Units','normalized',...
    'Position',[.45,.02,.1,0.08],'FontSize',14,...
    'Callback',{@Nextbutton_Callback});
    function Nextbutton_Callback(source,eventdata)    
        val=get(FRAPfileslistboxh,'Value');
        FileLocation=get(Filelocationh,'String');
        FileNames=get(FRAPfileslistboxh,'String');
        for j=1:length(val)
            fullfilepaths{j,1}=fullfile(FileLocation,FileNames{val(j)});
            basicinput{j,1}=FileNames{val(j)};
        end
        basicinput{1,2}=(get(FRAPmodelh,'Value'));
        basicinput{1,3}=(get(ROIh,'Value'));
        basicinput{1,4}=(get(normalizationh,'Value'));
        basicinput{1,5}=str2num(get(BGh,'String'));
        if isempty(basicinput{1,5})
            basicinput{1,5}=0;
            set(BGh,'String',0);
            warndlg('Background is incorrect, resetting to 0')
        end
        basicinput{1,6}=str2num(get(PBInumh,'String'));
        if isempty(basicinput{1,6}) | basicinput{1,6}<2
            basicinput{1,6}=2;
            set(PBInumh,'String',2);
            warndlg('Post bleach image number is incorrect, resetting to 2')
        end
        basicinput{1,7}=[0 0 0];
        basicinput{1,8}=str2num(get(PBNh,'String'));
        if isempty(basicinput{1,8}) | basicinput{1,8}<1
            basicinput{1,8}=1;
            set(PBNh,'String',1);
            warndlg('Number of images for normalization is incorrect, resetting to 1')
        end
        basicinput{1,9}=0;
        
        if basicinput{1,2}==1
            %% If Diffusion model
            prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
            dlg_title = 'Input';
            num_lines = 1;
            def = {'256 45','9'};
            answer = inputdlg(prompt,'ROI Input',num_lines,def);
            basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
            answer = questdlg('Do you want to calculate MF using an adjacent ROI?','MF calculation','No','Yes','Yes');
            switch answer
                case 'Yes'
                    basicinput{1,9} = 1;
                case 'No'
                    basicinput{1,9}=0;
            end
            
        try
            [data]=loadData_Diffusion(fullfilepaths,basicinput); %This performs image processing to obtain processed FRAP curves.
        catch errObj
%             if strcmp(errObj.identifier,'MATLAB:Java:GenericException');
%             errordlg('Bioformats does not support the filetype you selected.  Please select a Bioformats-supported image type.');     
            errordlg(getReport(errObj,'extended','hyperlinks','off'),'Error 1');        
            return
            %else
                %rethrow(errObj)
            %end
        end 
        
        
            assignin('base', 'data', data);
            assignin('base', 'basicinput', basicinput);
            assignin('base','FileLocation',FileLocation);
            Figure_GUI_Diffusion(data,basicinput) % This is the data analysis, visualization, and saving steps.
            
            
        elseif basicinput{1,2}==2
            %% If Reaction 1 model
            if basicinput{1,3}==1 % IF circular ROI
                prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'256 45','9'};
                answer = inputdlg(prompt,'ROI Input',num_lines,def);
                basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
                
            end
            [data]=loadData_Reaction(fullfilepaths,basicinput); %This performs image processing to obtain processed FRAP curves.
            assignin('base', 'data', data);
            assignin('base', 'basicinput', basicinput);
            assignin('base','FileLocation',FileLocation);
            Figure_GUI_Reaction(data,basicinput) % This is the data analysis, visualization, and saving steps.
            
            
        elseif basicinput{1,2}==3
            %% If Reaction 2 model
            if basicinput{1,3}==1 % IF circular ROI
                prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'256 45','9'};
                answer = inputdlg(prompt,'ROI Input',num_lines,def);
                basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
                
            end
            [data]=loadData_Reaction(fullfilepaths,basicinput); %This performs image processing to obtain processed FRAP curves.
            assignin('base', 'data', data);
            assignin('base', 'basicinput', basicinput);
            assignin('base','FileLocation',FileLocation);
            Figure_GUI_Reaction2(data,basicinput) % This is the data analysis, visualization, and saving steps.
            
        elseif basicinput{1,2}==4
            %% If FRAP/FRET model
            prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
            dlg_title = 'Input';
            num_lines = 1;
            def = {'256 45','9'};
            answer = inputdlg(prompt,'ROI Input',num_lines,def);
            basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
            answer = questdlg('Do you want to calculate MF using an adjacent ROI?','MF calculation','No','Yes','Yes');
            switch answer
                case 'Yes'
                    basicinput{1,9} = 1;
                case 'No'
                    basicinput{1,9}=0;
            end
            [data]=loadData_FRAP_FRET(fullfilepaths,basicinput);
            assignin('base', 'data', data);
            assignin('base', 'basicinput', basicinput);
            assignin('base','FileLocation',FileLocation);
            Figure_GUI_FRAP_FRET(data,basicinput)
            
        end
        
    end

Nextbuttonh = uicontrol('Parent',GUIfigureh,'Style','pushbutton',...
    'String','Preview','Units','normalized',...
    'Position',[.75,.02,.1,0.08],'FontSize',14,...
    'Callback',{@Previewbutton_Callback});
    function Previewbutton_Callback(source,eventdata)
        val=get(FRAPfileslistboxh,'Value');
        FileLocation=get(Filelocationh,'String');
        FileNames=get(FRAPfileslistboxh,'String');
        fullfilepath=fullfile(FileLocation,FileNames{val(1)});
        basicinput{1,1}=FileNames{val(1)};
        basicinput{1,2}=(get(FRAPmodelh,'Value'));
        basicinput{1,3}=(get(ROIh,'Value'));
        basicinput{1,6}=str2num(get(PBInumh,'String'));
        
        if basicinput{1,2}==1
            %% If Diffusion model
            prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
            dlg_title = 'Input';
            num_lines = 1;
            def = {'256 45','9'};
            answer = inputdlg(prompt,'ROI Input',num_lines,def);
            basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
            answer = questdlg('Do you want to calculate MF using an adjacent ROI?','MF calculation','No','Yes','Yes');
            switch answer
                case 'Yes'
                    basicinput{1,9} = 1;
                case 'No'
                    basicinput{1,9}=0;
            end
            imgdata=bfopen(fullfilepath);
            PreviewGUI_Diffusion(imgdata,basicinput);
            
        elseif basicinput{1,2}==2
            %% If Reaction 1 model
            imgdata=bfopen(fullfilepath);
            PreviewGUI_Reaction(imgdata,basicinput);
            
        elseif basicinput{1,2}==3
            %% If Reaction 2 model
            imgdata=bfopen(fullfilepath);
            PreviewGUI_Reaction(imgdata,basicinput);
            
        elseif basicinput{1,2}==4 % If FRAP/FRET model
            prompt = {'Pixel coordinates of ROI center','Pixel radius of ROI'};
            dlg_title = 'Input';
            num_lines = 1;
            def = {'256 45','9'};
            answer = inputdlg(prompt,'ROI Input',num_lines,def);
            basicinput{1,7}=[str2num(answer{1,1}),str2num(answer{2,1})];
            answer = questdlg('Do you want to calculate MF using an adjacent ROI?','MF calculation','No','Yes','Yes');
            switch answer
                case 'Yes'
                    basicinput{1,9} = 1;
                case 'No'
                    basicinput{1,9}=0;
            end
            imgdata=bfopen(fullfilepath);
            PreviewGUI_FRAP_FRET(imgdata,basicinput);
            
        end
        
        
    end

%% Initialize the GUI
movegui('center')
set(GUIfigureh,'Visible','on')

end
