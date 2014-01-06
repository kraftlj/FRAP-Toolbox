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

function [nucleusroimask]=ROIinitialization_NCtransport(img,basicinput)
% Inputs: img=image used for user defining ROIs; basicinput=basic user
% inputs from the main GUI.  Outputs: bleachroimask=mask for bleaching
% region; cellroimask=mask for cell region.

% User defined ROIs are required for the NCtransport model

h=msgbox({'Define the nuclear region using the mouse.  Right click and press create mask when finished.'});
movegui(h,'northwest');
figure
nucleusroimask = roipoly(imadjust(img,stretchlim(img,[0,1]))); % Create a mask using the user defined ROI
text('Position',[.9,.9],'String',{'Define the nuclear region using the mouse.  Right click and press create mask when finished.'});
h2=gcf;
close(h);
close(h2);


end