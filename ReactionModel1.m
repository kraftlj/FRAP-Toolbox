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

function [data, avg]=ReactionModel1(data, basicinput, usrinputs, val)
% Inputs:
% data - this is the data output from PhotoDecay_Reaction
% basicinput - this is basic user input from the Main_GUI
% usrinputs - this is basic user input from Figure_GUI_Reaction
% val - the datasets that were selected by the user for plotting

options=optimset('lsqcurvefit');
options.Display='off';
%% Fit the averaged data;
if usrinputs{7,1}==2
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
    end
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(1).time(basicinput{1,6}+usrinputs{5,1}-1);
    avg.t=t+t(usrinputs{5,1});
    %----------------------------------------------------------------------
    
    wf=f./(t+sum(f)); % Weight the beginning of the FRAP curve more
    switch [usrinputs{1,4},usrinputs{2,4},usrinputs{3,4}]
        case 'AdjustableAdjustableAdjustable'
            frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3}],options);
            avg.a=p(1);
            avg.b=p(2);
            avg.c=p(3);
        case 'FixedAdjustableAdjustable'
            avg.a=usrinputs{1,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{2,2},usrinputs{3,2}],[usrinputs{2,3},usrinputs{3,3}],options);
            avg.b=p(1);
            avg.c=p(2);
        case 'AdjustableFixedAdjustable'
            avg.b=usrinputs{2,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{3,3}],options);
            avg.a=p(1);
            avg.c=p(2);
        case 'AdjustableAdjustableFixed'
            avg.c=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1}],t,wf,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
            avg.a=p(1);
            avg.b=p(2);
        case 'FixedFixedAdjustable'
            avg.a=usrinputs{1,1};
            avg.b=usrinputs{2,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
            avg.c=p(1);
        case 'FixedAdjustableFixed'
            avg.a=usrinputs{1,1};
            avg.c=usrinputs{3,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
            avg.b=p(1);
        case 'AdjustableFixedFixed'
            avg.b=usrinputs{2,1};
            avg.c=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
            avg.a=p(1);
        case 'FixedFixedFixed'
            avg.a=usrinputs{1,1};
            avg.b=usrinputs{2,1};
            avg.c=usrinputs{3,1};
    end
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
    frapfit=frapfun([avg.a,avg.b,avg.c],t).*(t+sum(f)); % unweight the optimized fit
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([avg.a,avg.b,avg.c],t); % weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of square errors
    %----------------------------------------------------------------------
    
    %% Fit the individual FRAP datasets
else
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        t=data(val(index1)).time(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(val(index1)).time(basicinput{1,6}+usrinputs{5,1}-1);
        f=data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
        data(val(index1)).t=t+t(usrinputs{5,1});
        wf=f./(t+sum(f)); % weight the beginning of the FRAP curve more
        switch [usrinputs{1,4},usrinputs{2,4},usrinputs{3,4}]
            case 'AdjustableAdjustableAdjustable'
                frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3}],options);
                data(val(index1)).a=p(1);
                data(val(index1)).b=p(2);
                data(val(index1)).c=p(3);
            case 'FixedAdjustableAdjustable'
                data(val(index1)).a=usrinputs{1,1};
                frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{2,2},usrinputs{3,2}],[usrinputs{2,3},usrinputs{3,3}],options);
                data(val(index1)).b=p(1);
                data(val(index1)).c=p(2);
            case 'AdjustableFixedAdjustable'
                data(val(index1)).b=usrinputs{2,1};
                frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{3,3}],options);
                data(val(index1)).a=p(1);
                data(val(index1)).c=p(2);
            case 'AdjustableAdjustableFixed'
                data(val(index1)).c=usrinputs{3,1};
                frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1}],t,wf,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
                data(val(index1)).a=p(1);
                data(val(index1)).b=p(2);
            case 'FixedFixedAdjustable'
                data(val(index1)).a=usrinputs{1,1};
                data(val(index1)).b=usrinputs{2,1};
                frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
                data(val(index1)).c=p(1);
            case 'FixedAdjustableFixed'
                data(val(index1)).a=usrinputs{1,1};
                data(val(index1)).c=usrinputs{3,1};
                frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
                data(val(index1)).b=p(1);
            case 'AdjustableFixedFixed'
                data(val(index1)).b=usrinputs{2,1};
                data(val(index1)).c=usrinputs{3,1};
                frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)) ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
                data(val(index1)).a=p(1);
            case 'FixedFixedFixed'
                data(val(index1)).a=usrinputs{1,1};
                data(val(index1)).b=usrinputs{2,1};
                data(val(index1)).c=usrinputs{3,1};
        end
        frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
        frapfit=frapfun([data(val(index1)).a,data(val(index1)).b,data(val(index1)).c],t).*(t+sum(f)); %unweight the optimized FRAP fit
        data(val(index1)).frapfit=frapfit;
        data(val(index1)).frapres=wf-frapfun([data(val(index1)).a,data(val(index1)).b,data(val(index1)).c],t); % weighted residuals
        data(val(index1)).SS=sum(data(val(index1)).frapres.^2); % Sum of squared errors
    end
    
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(1).time(basicinput{1,6}+usrinputs{5,1}-1);
    wf=f./(t+sum(f));
    avg.t=t+t(usrinputs{5,1});
    a=mean([data.a]);
    avg.a=a;
    b=mean([data.b]);
    avg.b=b;
    c=mean([data.c]);
    avg.c=c;
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
    frapfit=frapfun([a,b,c],t).*(t+sum(f)); % unweighted optimized fit using the averaged parameters
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([a,b,c],t); % Weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of squared errors
end

end
