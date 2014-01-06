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

function [data]=loadData_NCtransport(fullfilepaths,basicinput)
% Inputs: basicinput=basic user input from main GUI; fullfilepaths=cell
% array of the full paths to user selected files in the main GUI.
% Outputs:
% Function initialization--------------------------------------------------
h=waitbar(0/length(fullfilepaths),['Loading Data (Initializing)']);
for index1=1:length(fullfilepaths) % Loop through each file one by one
    waitbar((index1-1)/length(fullfilepaths),h,['Loading Data: ', num2str((index1-1)/length(fullfilepaths)*100),'%'])
    
    [pathstr, name, ext] = fileparts(fullfilepaths{index1}); % define the file parts
    
    imgdata=bfopen(fullfilepaths{index1}); % Load the image file using open Bio-Formats
    
%     metadata=imgdata{1,2}; % Load the meta data information into variable metadata
    omeMeta=imgdata{1,4};
    img=imgdata{1,1}{1,1}; % Load the first image plane in the stack for use in user defining ROI
    % End of Function initialization-------------------------------------------
    
    % ROI initialization-------------------------------------------------------
    [nucleusroimask]=ROIinitialization_NCtransport2(img,basicinput);
    % End of ROI intialization-------------------------------------------------
    
    % Record FRAP curve information--------------------------------------------
    [data(index1).time data(index1).nucleus]=FRAPcurve_NCtransport2(imgdata,basicinput,omeMeta,nucleusroimask,name);
    % End of Record FRAP curve information-------------------------------------
    
    % Normalize the FRAP curve-------------------------------------------------
    [data(index1).normfrap]=NormalizeFRAP_NCtransport2(data(index1).nucleus,basicinput);
    % End of Normalize the FRAP curve------------------------------------------
    
end
close(h)
end
