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

function [data, avg]=ReactionModel2(data, basicinput, usrinputs, val)
% Inputs:
% data - this is the data output from PhotoDecay_Reaction
% basicinput - this is basic user input from the Main_GUI
% usrinputs - this is basic user input from Figure_GUI_Reaction
% val - the datasets that were selected by the user for plotting

options=optimset('lsqcurvefit');
options.Display='off';
%% Fit the averaged data;
if usrinputs{8,1}==2 
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
    end
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2})-data(1).time(basicinput{1,6}+usrinputs{6,1}-1);
    avg.t=t+t(usrinputs{6,1});
    %----------------------------------------------------------------------
    
    wf=f./(t+sum(f)); % Weight the beginning of the FRAP curve more
    switch [usrinputs{1,4},usrinputs{2,4},usrinputs{3,4},usrinputs{4,4}]
        case 'AdjustableAdjustableAdjustableAdjustable'
            frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-p(4)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3},usrinputs{4,3}],options);
            avg.Finf=p(1);
            avg.C1eq=p(2);
            avg.koff1=p(3);
            avg.koff2=p(4);
        case 'FixedAdjustableAdjustableAdjustable'
            avg.Finf=usrinputs{1,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{2,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{2,3},usrinputs{3,3},usrinputs{4,3}],options);
            avg.C1eq=p(1);
            avg.koff1=p(2);
            avg.koff2=p(3);
        case 'AdjustableFixedAdjustableAdjustable'
            avg.C1eq=usrinputs{2,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{3,3},usrinputs{4,3}],options);
            avg.Finf=p(1);
            avg.koff1=p(2);
            avg.koff2=p(3);
        case 'AdjustableAdjustableFixedAdjustable'
            avg.koff1=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-p(2))*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{4,3}],options);
            avg.Finf=p(1);
            avg.C1eq=p(2);
            avg.koff2=p(3);
        case 'AdjustableAdjustableAdjustableFixed'
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3}],options);
            avg.Finf=p(1);
            avg.C1eq=p(2);
            avg.koff1=p(3);
        case 'FixedFixedAdjustableAdjustable'
            avg.Finf=usrinputs{1,1};
            avg.C1eq=usrinputs{2,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{3,2},usrinputs{4,2}],[usrinputs{3,3},usrinputs{4,3}],options);
            avg.koff1=p(1);
            avg.koff2=p(2); 
        case 'FixedAdjustableFixedAdjustable'
            avg.Finf=usrinputs{1,1};
            avg.koff1=usrinputs{3,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{4,1}],t,wf,[usrinputs{2,2},usrinputs{4,2}],[usrinputs{2,3},usrinputs{4,3}],options);
            avg.C1eq=p(1);
            avg.koff2=p(2);
        case 'FixedAdjustableAdjustableFixed'
            avg.Finf=usrinputs{1,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{2,2},usrinputs{3,2}],[usrinputs{2,3},usrinputs{3,3}],options);
            avg.C1eq=p(1);
            avg.koff1=p(2);
        case 'AdjustableFixedFixedAdjustable'
            avg.C1eq=usrinputs{2,1};
            avg.koff1=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{4,3}],options);
            avg.Finf=p(1);
            avg.koff2=p(2);
        case 'AdjustableFixedAdjustableFixed'
            avg.C1eq=usrinputs{2,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{3,3}],options);
            avg.Finf=p(1);
            avg.koff1=p(2);    
        case 'AdjustableAdjustableFixedFixed'
            avg.koff1=usrinputs{3,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-p(2))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1}],t,wf,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
            avg.Finf=p(1);
            avg.C1eq=p(2);             
        case 'FixedFixedFixedAdjustable'
            avg.Finf=usrinputs{1,1};
            avg.C1eq=usrinputs{2,1};
            avg.koff1=usrinputs{3,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-p(1)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{4,1}],t,wf,[usrinputs{4,2}],[usrinputs{4,3}],options);
            avg.koff2=p(1);  
        case 'FixedFixedAdjustableFixed'
            avg.Finf=usrinputs{1,1};
            avg.C1eq=usrinputs{2,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
            avg.koff1=p(1); 
        case 'FixedAdjustableFixedFixed'
            avg.Finf=usrinputs{1,1};
            avg.koff1=usrinputs{3,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
            avg.C1eq=p(1);  
        case 'AdjustableFixedFixedFixed'
            avg.koff1=usrinputs{3,1};
            avg.C1eq=usrinputs{2,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
            avg.Finf=p(1);             
        case 'FixedFixedFixedFixed'
            avg.Finf=usrinputs{1,1};
            avg.koff1=usrinputs{3,1};
            avg.C1eq=usrinputs{2,1};
            avg.koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
               
    end
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-p(4)*t)) ./(t+sum(f));
    frapfit=frapfun([avg.Finf,avg.C1eq,avg.koff1,avg.koff2],t).*(t+sum(f)); % unweight the optimized fit
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([avg.Finf,avg.C1eq,avg.koff1,avg.koff2],t); % weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of square errors
    %----------------------------------------------------------------------

    %% Fit the individual FRAP datasets    
else
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        t=data(val(index1)).time(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2})-data(val(index1)).time(basicinput{1,6}+usrinputs{6,1}-1);
        f=data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2});
        data(val(index1)).t=t+t(usrinputs{6,1});
        wf=f./(t+sum(f)); % weight the beginning of the FRAP curve more
        switch [usrinputs{1,4},usrinputs{2,4},usrinputs{3,4},usrinputs{4,4}]
        case 'AdjustableAdjustableAdjustableAdjustable'
            frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-p(4)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3},usrinputs{4,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).C1eq=p(2);
            data(val(index1)).koff1=p(3);
            data(val(index1)).koff2=p(4);
        case 'FixedAdjustableAdjustableAdjustable'
            data(val(index1)).Finf=usrinputs{1,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{2,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{2,3},usrinputs{3,3},usrinputs{4,3}],options);
            data(val(index1)).C1eq=p(1);
            data(val(index1)).koff1=p(2);
            data(val(index1)).koff2=p(3);
        case 'AdjustableFixedAdjustableAdjustable'
            data(val(index1)).C1eq=usrinputs{2,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{3,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{3,3},usrinputs{4,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).koff1=p(2);
            data(val(index1)).koff2=p(3);
        case 'AdjustableAdjustableFixedAdjustable'
            data(val(index1)).koff1=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-p(2))*exp(-p(3)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{4,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).C1eq=p(2);
            data(val(index1)).koff2=p(3);
        case 'AdjustableAdjustableAdjustableFixed'
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).C1eq=p(2);
            data(val(index1)).koff1=p(3);
        case 'FixedFixedAdjustableAdjustable'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).C1eq=usrinputs{2,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{3,2},usrinputs{4,2}],[usrinputs{3,3},usrinputs{4,3}],options);
            data(val(index1)).koff1=p(1);
            data(val(index1)).koff2=p(2); 
        case 'FixedAdjustableFixedAdjustable'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).koff1=usrinputs{3,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{4,1}],t,wf,[usrinputs{2,2},usrinputs{4,2}],[usrinputs{2,3},usrinputs{4,3}],options);
            data(val(index1)).C1eq=p(1);
            data(val(index1)).koff2=p(2);
        case 'FixedAdjustableAdjustableFixed'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-p(2)*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{2,2},usrinputs{3,2}],[usrinputs{2,3},usrinputs{3,3}],options);
            data(val(index1)).C1eq=p(1);
            data(val(index1)).koff1=p(2);
        case 'AdjustableFixedFixedAdjustable'
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff1=usrinputs{3,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-p(2)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{4,1}],t,wf,[usrinputs{1,2},usrinputs{4,2}],[usrinputs{1,3},usrinputs{4,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).koff2=p(2);
        case 'AdjustableFixedAdjustableFixed'
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-p(2)*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{3,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).koff1=p(2);    
        case 'AdjustableAdjustableFixedFixed'
            data(val(index1)).koff1=usrinputs{3,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-p(2)*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-p(2))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1}],t,wf,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
            data(val(index1)).Finf=p(1);
            data(val(index1)).C1eq=p(2);             
        case 'FixedFixedFixedAdjustable'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff1=usrinputs{3,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-p(1)*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{4,1}],t,wf,[usrinputs{4,2}],[usrinputs{4,3}],options);
            data(val(index1)).koff2=p(1);  
        case 'FixedFixedAdjustableFixed'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-p(1)*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
            data(val(index1)).koff1=p(1); 
        case 'FixedAdjustableFixedFixed'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).koff1=usrinputs{3,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-p(1)*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-p(1))*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
            data(val(index1)).C1eq=p(1);  
        case 'AdjustableFixedFixedFixed'
            data(val(index1)).koff1=usrinputs{3,1};
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (p(1)-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(p(1)-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
            data(val(index1)).Finf=p(1);             
        case 'FixedFixedFixedFixed'
            data(val(index1)).Finf=usrinputs{1,1};
            data(val(index1)).koff1=usrinputs{3,1};
            data(val(index1)).C1eq=usrinputs{2,1};
            data(val(index1)).koff2=usrinputs{4,1};
            frapfun=@(p,t) (usrinputs{1,1}-usrinputs{2,1}*exp(-usrinputs{3,1}*t)-(usrinputs{1,1}-f(1)-usrinputs{2,1})*exp(-usrinputs{4,1}*t)) ./(t+sum(f));
               
    end
        frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-p(4)*t)) ./(t+sum(f));
        frapfit=frapfun([data(val(index1)).Finf,data(val(index1)).C1eq,data(val(index1)).koff1,data(val(index1)).koff2],t).*(t+sum(f)); %unweight the optimized FRAP fit
        data(val(index1)).frapfit=frapfit;
        data(val(index1)).frapres=wf-frapfun([data(val(index1)).Finf,data(val(index1)).C1eq,data(val(index1)).koff1,data(val(index1)).koff2],t); % weighted residuals
        data(val(index1)).SS=sum(data(val(index1)).frapres.^2); % Sum of squared errors
    end
    
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{6,1}-1:usrinputs{6,2})-data(1).time(basicinput{1,6}+usrinputs{6,1}-1);
    wf=f./(t+sum(f));
    avg.t=t+t(usrinputs{6,1});
    koff1=mean([data.koff1]);
    avg.koff1=koff1;
    Finf=mean([data.Finf]);
    avg.Finf=Finf;
    C1eq=mean([data.C1eq]);
    avg.C1eq=C1eq;
    koff2=mean([data.koff2]);
    avg.koff2=koff2;
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)-(p(1)-f(1)-p(2))*exp(-p(4)*t)) ./(t+sum(f));
    frapfit=frapfun([Finf,C1eq,koff1,koff2],t).*(t+sum(f)); % unweighted optimized fit using the averaged parameters
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([Finf,C1eq,koff1,koff2],t); % Weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of squared errors
end

end
