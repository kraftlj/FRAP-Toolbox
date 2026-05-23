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

function [normfrap]=NormalizeFRAP_Reaction(frap,cell,basicinput)
% Inputs:
% frap - mean intensity inside bleach region.
% cell - mean intensity inside cell region.
% basicinput - basic user inputs from the main GUI.

% Outputs:
% normfrap - normalized FRAP curve.

if basicinput{1,4}==2 % IF normalizing by whole cell intensity
    normfrap=frap./cell; % The cell normalization corrects for small movement in z-plane and photofading
    normfrap=normfrap./mean(normfrap(basicinput{1,6}-1-basicinput{1,8}+1:basicinput{1,6}-1)); %The FRAP data is normalized to 1.
        
else % IF only normalizing by the pre-bleach steady-state intensity
    normfrap=frap./mean(frap(basicinput{1,6}-1-basicinput{1,8}+1:basicinput{1,6}-1)); %The FRAP data is normalized to 1.
    
    
end

end