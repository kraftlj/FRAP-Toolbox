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

function PreviewGUI(imgdata,basicinput)
GUIfigureh = figure('Position',[200 200 1000 400],...
    'MenuBar','none','NumberTitle','off',...
    'Name','Preview of an Image Stack','Visible','off');

img=imgdata{1,1}{1,1};

if basicinput{1,3}==1 % IF circular ROI
    
    x0=basicinput{1,7}(1); % X coordinates for the center of a circular bleaching ROI
    y0=basicinput{1,7}(2); % Y coordinates for the center of a circular bleaching ROI
    R0=basicinput{1,7}(3); % Radius of circular bleaching ROI
    t = 0:pi/20:2*pi;
    xi = R0*cos(t)+x0; % X coordinates for the circle
    yi = R0*sin(t)+y0; % Y coordinates for the circle
    bleachroimask = poly2mask(xi,yi, size(img,1),size(img,2));
    if basicinput{1,9}==1;
        xi = R0*cos(t)+x0+R0*2.5; % X coordinates for the circle
        yi = R0*sin(t)+y0; % Y coordinates for the circle
        adjacentroimask = poly2mask(xi,yi, size(img,1),size(img,2));
        BW2 = bwperim(bleachroimask);
        SE = strel('disk',1);
        BW2=imdilate(BW2,SE);
        img(BW2)=max(img(:));
        BW3 = bwperim(adjacentroimask);
        BW3=imdilate(BW3,SE);
        img(BW3)=max(img(:));
    else
        BW2 = bwperim(bleachroimask);
        SE = strel('disk',1);
        BW2=imdilate(BW2,SE);
        img(BW2)=max(img(:));
    end
    
end
imagesc(img)
axis image
colormap(gray)
title(basicinput{1,1});

imgnumberh=uicontrol('Parent',GUIfigureh,'Style','text','String','1',...
    'Units','normalized','Position',[.05,.02,.05,.1],'FontSize',14,'BackgroundColor',get(GUIfigureh,'color'));

maxd=length(imgdata{1,1});
sliderh=uicontrol('Parent',GUIfigureh,'Style','slider','Min',1,'Max',maxd,'SliderStep',[1/(maxd-1),(1/(maxd-1))*10],'Value',1,...,
    'Units','normalized','Position',[.1,.02,.8,.1],'FontSize',14,'Callback',{@Slider_Callback});
    function Slider_Callback(source,eventdata)
        val=round(get(sliderh,'Value'));
        assignin('base','val',val);
        img=imgdata{1,1}{val,1};
        if basicinput{1,3}==1 % IF circular ROI
            img(BW2)=max(img(:));
            if basicinput{1,9}==1
                img(BW3)=max(img(:));
            end
        end
        imagesc(img);
        axis image
        colormap(gray)
        title(basicinput{1,1});
        set(imgnumberh,'String',num2str(val));
        
        
    end

%--------------------------------------------------------------------------
% Move the GUI to the center of the screen.
movegui(GUIfigureh,'center')
% Make the GUI visible.
set(GUIfigureh,'Visible','on');

%--------------------------------------------------------------------------

end