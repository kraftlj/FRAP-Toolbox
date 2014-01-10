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

function [data, avg]=DiffusionModel(data, basicinput, usrinputs, val)
% Inputs:
% data - this is the data output from PhotoDecay.m
% basicinput - these are basic user inputs from the Main_GUI
% usrinputs - these are basic user inputs from Figure_GUI_Diffusion
% val - these are the FRAP datasets selected by the user for analysis.

% Outputs:
% data - this is the input data amended with the optimized parameters and
% the optimized FRAP curves and initial conditions.
% avg - this is the averaged frap data, optimized parameters, frap curves,
% and initial conditions.

%%
options=optimset('lsqcurvefit');
options.Display='off';
if usrinputs{10,1}==2  % Fit the averaged data
    avgpbp=[];
    avgr=[];
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        avgadjacent(index1,:)=data(val(index1)).adjacent;
        avgpbp(end+1:end+length(data(val(index1)).pbp))=data(val(index1)).pbp;
        avgr(end+1:end+length(data(val(index1)).r))=data(val(index1)).r;
    end
    adjacent=mean(avgadjacent,1);
    f=mean(avgf,1);
    avg.f=f;
    avg.correctMF=1-(mean(adjacent(usrinputs{9,1}:usrinputs{9,2}))-mean(f(usrinputs{9,1}:usrinputs{9,2})));
    f=f(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2});
    A=[avgr',avgpbp']; % Average the data at equivalent radial distances from the center of the bleach ROI.
    [u,~,id2] = unique(A(:,1),'rows');
    B = [u,accumarray(id2,A(:,2))./accumarray(id2,1)];
    r=B(:,1);
    avg.r=r;
    pbp=B(:,2);
    avg.pbp=pbp;
    pbp=pbp(usrinputs{6,1}:usrinputs{6,2});
    r=r(usrinputs{6,1}:usrinputs{6,2});
    voxelSizeX=data(1).voxelSizeX;
    rn=basicinput{1,7}(3).*voxelSizeX;
    avg.rfit=r;
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2})-data(1).time(basicinput{1,6});
    avg.t=t;
    
    %% Find the optimized initial conditions
    switch [usrinputs{1,4},usrinputs{2,4}]
        case 'AdjustableAdjustable'
            pbpfun=@(p,r) exp(-p(1)*exp(-2*r.^2/p(2)^2));
            p=lsqcurvefit(pbpfun,[usrinputs{1,1},usrinputs{2,1}],r,pbp,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
            k=p(1);
            re=p(2);
            frapfun=@(k,t) KangFRAP(t,re,rn,1,k);
            k=lsqcurvefit(frapfun,[1],0,f(1),[0],[Inf],options);
            avg.k=k;
            pbpfun=@(re,r) exp(-k*exp(-2*r.^2/re^2));
            re=lsqcurvefit(pbpfun,[1],r,pbp,[0],[Inf],options);
            avg.re=re;
        case 'FixedAdjustable'
            k=usrinputs{1,1};
            pbpfun=@(p,r) exp(-k*exp(-2*r.^2/p(1)^2));
            p=lsqcurvefit(pbpfun,[usrinputs{2,1}],r,pbp,[usrinputs{2,2}],[usrinputs{2,3}],options);
            re=p(1);
            avg.k=k;
            avg.re=re;
        case 'AdjustableFixed'
            re=usrinputs{2,1};
            frapfun=@(k,t) KangFRAP(t,re,rn,1,k);
            k=lsqcurvefit(frapfun,[1],0,f(1),[0],[Inf],options);
            avg.k=k;
            avg.re=re;
        case 'FixedFixed'
            k=usrinputs{1,1};
            avg.k=k;
            re=usrinputs{2,1};
            avg.re=re;
    end
    
    pbpfun=@(p,r) exp(-p(1)*exp(-2*r.^2/p(2)^2));
    pbpfit=pbpfun([k,re],r);
    avg.pbpfit=pbpfit;
    avg.pbpres=pbp-pbpfit;
    
    %% Find the optimized FRAP curves
    wf=f./(t+sum(f)); % Weight the beginning of the FRAP data more heavily
    switch [usrinputs{3,4},usrinputs{4,4}]
        case 'AdjustableAdjustable'
            frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*p(2)+(1-p(2))*f(1))./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{3,2},usrinputs{4,2}],[usrinputs{3,3},usrinputs{4,3}],options);
            avg.D=p(1);
            avg.MF=p(2);
        case 'FixedAdjustable'
            avg.D=usrinputs{3,1};
            frapfun=@(p,t) (KangFRAP(t,re,rn,avg.D,k).*p(1)+(1-p(1))*f(1))./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{4,1}],t,wf,[usrinputs{4,2}],[usrinputs{4,3}],options);
            avg.MF=p(1);
        case 'AdjustableFixed'
            avg.MF=usrinputs{4,1};
            frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*avg.MF+(1-avg.MF)*f(1))./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
            avg.D=p(1);
        case 'FixedFixed'
            avg.D=usrinputs{3,1};
            avg.MF=p(1);
    end
    
    frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*p(2)+(1-p(2))*f(1))./(t+sum(f));
    frapfit=frapfun([avg.D,avg.MF],t).*(t+sum(f));
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([avg.D,avg.MF],t);
    avg.SS=sum([avg.frapres].^2);
    
    %--------------------------------------------------------------------------
    
    %% Fit the each individual FRAP dataset independently.
else
    avgpbp=[];
    avgr=[];
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        avgadjacent(index1,:)=data(val(index1)).adjacent;
        avgpbp(end+1:end+length(data(val(index1)).pbp))=data(val(index1)).pbp;
        avgr(end+1:end+length(data(val(index1)).r))=data(val(index1)).r;
        data(val(index1)).correctMF=1-(mean(data(val(index1)).adjacent(usrinputs{9,1}:usrinputs{9,2}))-mean(data(val(index1)).normfrap(usrinputs{9,1}:usrinputs{9,2})));
        voxelSizeX=data.voxelSizeX;
        rn=basicinput{1,7}(3).*voxelSizeX;
        r=data(val(index1)).r(usrinputs{6,1}:usrinputs{6,2});
        data(val(index1)).rfit=r;
        pbp=data(val(index1)).pbp(usrinputs{6,1}:usrinputs{6,2});
        t=data(val(index1)).time(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2})-data(val(index1)).time(basicinput{1,6});
        f=data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2});
        data(val(index1)).t=t;
        %% Find the optimized initial conditions
        switch [usrinputs{1,4},usrinputs{2,4}]
            case 'AdjustableAdjustable'
                pbpfun=@(p,r) exp(-p(1)*exp(-2*r.^2/p(2)^2));
                p=lsqcurvefit(pbpfun,[usrinputs{1,1},usrinputs{2,1}],r,pbp,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
                k=p(1);
                re=p(2);
                frapfun=@(k,t) KangFRAP(t,re,rn,1,k);
                k=lsqcurvefit(frapfun,[1],0,f(1),[0],[Inf],options);
                data(val(index1)).k=k;
                pbpfun=@(re,r) exp(-k*exp(-2*r.^2/re^2));
                re=lsqcurvefit(pbpfun,[1],r,pbp,[0],[Inf],options);
                data(val(index1)).re=re;
            case 'FixedAdjustable'
                k=usrinputs{1,1};
                pbpfun=@(p,r) exp(-k*exp(-2*r.^2/p(1)^2));
                p=lsqcurvefit(pbpfun,[usrinputs{2,1}],r,pbp,[usrinputs{2,2}],[usrinputs{2,3}],options);
                re=p(1);
                data(val(index1)).k=k;
                data(val(index1)).re=re;
            case 'AdjustableFixed'
                re=usrinputs{2,1};
                frapfun=@(k,t) KangFRAP(t,re,rn,1,k);
                k=lsqcurvefit(frapfun,[1],0,f(1),[0],[Inf],options);
                data(val(index1)).k=k;
                data(val(index1)).re=re;
            case 'FixedFixed'
                k=usrinputs{1,1};
                data(val(index1)).k=k;
                re=usrinputs{2,1};
                data(val(index1)).re=re;
        end
        
        pbpfun=@(p,r) exp(-p(1)*exp(-2*r.^2/p(2)^2));
        pbpfit=pbpfun([k,re],r);
        data(val(index1)).pbpfit=pbpfit;
        data(val(index1)).pbpres=pbp-pbpfit;
        
        %% Find the optimized FRAP curves
        wf=f./(t+sum(f)); % Weight the beginning of the FRAP datasets more heavily
        switch [usrinputs{3,4},usrinputs{4,4}]
            case 'AdjustableAdjustable'
                frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*p(2)+(1-p(2))*f(1))./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{3,1},usrinputs{4,1}],t,wf,[usrinputs{3,2},usrinputs{4,2}],[usrinputs{3,3},usrinputs{4,3}],options);
                data(val(index1)).D=p(1);
                data(val(index1)).MF=p(2);
            case 'FixedAdjustable'
                data(val(index1)).D=usrinputs{3,1};
                frapfun=@(p,t) (KangFRAP(t,re,rn,data(val(index1)).D,k).*p(1)+(1-p(1))*f(1))./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{4,1}],t,wf,[usrinputs{4,2}],[usrinputs{4,3}],options);
                data(val(index1)).MF=p(1);
            case 'AdjustableFixed'
                data(val(index1)).MF=usrinputs{4,1};
                frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*data(val(index1)).MF+(1-data(val(index1)).MF)*f(1))./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{3,1}],t,wf,[usrinputs{3,2}],[usrinputs{3,3}],options);
                data(val(index1)).D=p(1);
            case 'FixedFixed'
                data(val(index1)).D=usrinputs{3,1};
                data(val(index1)).MF=p(1);
        end
        
        frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*p(2)+(1-p(2))*f(1))./(t+sum(f));
        frapfit=frapfun([data(val(index1)).D,data(val(index1)).MF],t).*(t+sum(f));
        data(val(index1)).frapfit=frapfit;
        data(val(index1)).frapres=wf-frapfun([data(val(index1)).D,data(val(index1)).MF],t);
        data(val(index1)).SS=sum([data(val(index1)).frapres].^2);
    end
    % Average the parameters from the individual fits.
    f=mean(avgf,1);
    avg.f=f;
    avg.correctMF=mean([data.correctMF]);
    f=f(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2});
    A=[avgr',avgpbp'];
    [u,~,id2] = unique(A(:,1),'rows');
    B = [u,accumarray(id2,A(:,2))./accumarray(id2,1)];
    r=B(:,1);
    avg.r=r;
    pbp=B(:,2);
    avg.pbp=pbp;
    rfit=r(usrinputs{6,1}:usrinputs{6,2});
    avg.rfit=rfit;
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{7,1}-1:usrinputs{7,2})-data(1).time(basicinput{1,6});
    wf=f./(t+sum(f));
    avg.t=t;
    k=mean([data.k]);
    avg.k=k;
    re=mean([data.re]);
    avg.re=re;
    D=mean([data.D]);
    avg.D=D;
    MF=mean([data.MF]);
    avg.MF=MF;
    pbpfit=pbpfun([k,re],rfit);
    avg.pbpfit=pbpfit;
    avg.pbpres=pbp(usrinputs{6,1}:usrinputs{6,2})-pbpfit;
    frapfun=@(p,t) (KangFRAP(t,re,rn,p(1),k).*p(2)+(1-p(2))*f(1))./(t+sum(f));
    frapfit=frapfun([D,MF],t).*(t+sum(f));
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([D,MF],t);
    avg.SS=sum([avg.frapres].^2);
end

%% This is the theoretical diffusion FRAP function.

    function y=KangFRAP(t,re,rn,D,K)
        % Inputs:
        % t - the time information for the FRAP curves
        % re - the effective radius from the initial conditions
        % rn - the nominal radius of the user defined bleach region
        % D - the diffusion coefficient
        % K - the bleach depth from the initial conditions
        
        % Outputs:
        % y - the theoretical diffusion FRAP curve
        
        % Refer to Kang et al. 2008 for a discussion of this equation.
        % Partial summation from m=0 to m=10
        m=10;
        a=((-K).^(0:m)'*re^2) ./factorial(0:m)'*ones(1,length(t));
        b=re^2+((0:m)'*(8*D*t+rn^2));
        c=a./b;
        y=sum(c,1);
        
    end
end
