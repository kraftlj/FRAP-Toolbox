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

function [bleachroimask cellroimask]=ROIinitialization_Reaction(img,basicinput)
% Inputs:
% img - image used for user defining ROIs.
% basicinput - basic user inputs from the main GUI.

% Outputs:
% bleachroimask - mask for bleaching region.
% cellroimask - mask for cell region.

%%
if basicinput{1,3}==1 % If circular ROI
    x0=basicinput{1,7}(1); % X coordinates for the center of a circular bleaching ROI
    y0=basicinput{1,7}(2); % Y coordinates for the center of a circular bleaching ROI
    R0=basicinput{1,7}(3); % Radius of circular bleaching ROI
    [columnsInImage rowsInImage] = meshgrid(1:size(img,2), 1:size(img,1)); % Generate a meshgrid for the mask in the next steps
    bleachroimask = (rowsInImage - y0).^2 + (columnsInImage - x0).^2 <= R0.^2; % Create a mask using the circular coordinates
else % If using a user defined bleach region
    h=msgbox({'Define the bleach region using the mouse.  Right click and press create mask when finished.'});
    movegui(h,'northwest');
    figure
    bleachroimask = roipoly(imadjust(img,stretchlim(img,[0,1]))); % Create a mask using the user defined ROI
    h2=gcf;
    close(h);
    close(h2);
end

if basicinput{1,4}==2 % IF normalizing by whole cell intensity
    h=msgbox({'Define the cell region using the mouse.  Right click and press create mask when finished.'});
    movegui(h,'northwest');
    figure
    cellroimask = roipoly(imadjust(img,stretchlim(img,[0,1]))); % Create a mask using the user defined ROI
    h2=gcf;
    close(h);
    close(h2);
else
    cellroimask = NaN;
end

end