clear; close all;
load Phase.mat
currentgenpath= pwd;
subjectall = dir([currentgenpath,'\','S*.mat']);
% NFFT=2^nextpow2(1250);
Fs=250;
NFFT=10*Fs;
f=250/2.*linspace(0,1,NFFT/2+1);
freSeq=[8 12 11.2 15.2 10.4 14.4 9.6 13.6 8.8 12.8];
channel=62;
phaLoc=find(phases==0);

for fold_num=1:length(subjectall)
    tic
    subject=subjectall(fold_num).name;
    load ([currentgenpath,'\',subject])
    % calculate phase
    ytmp=data(channel,126:1375,phaLoc,:);
    Y=cell(10,6);phases=cell(10,6);
    Phase=[];SNR=[];
    for trial=1:10
        stif=freSeq(trial);
        phaCom=[];SNRCom=[];
        for block=1:6
            y=squeeze(ytmp(1,:,trial,block));
            [~,freLocInY]=min(abs(f-stif));[~,freLocInY_1]=min(abs(f-(stif-1)));[~,freLocInY_01]=min(abs(f-(stif+1)));
            % get SNR & phases for all fft frequencies
            [ frequency,Y{trial,block}] = fft_plot(y,Fs);
            phases{trial,block}=radtodeg(atan2(imag(Y{trial,block}),real(Y{trial,block})))-180; % get phase for all frequencies in fft
            phase(trial,block)=phases{trial,block}(freLocInY);
            phaComtmp=phase(trial,block);
            phaCom=[phaCom,phaComtmp];
            Y_amp=2*abs(Y{trial,block}(1:NFFT/2+1)); %amplitude after fft for SNR
            averIn1= mean(Y_amp(freLocInY_1:freLocInY_01)); %average amplitude of stimulus+ -1Hz
            SNR_indextmp= Y_amp(freLocInY)/averIn1;
            SNRCom=[SNRCom,SNR_indextmp];
        end
        Phase=[Phase;phaCom];
        SNR=[SNR;SNRCom];
    end
    
    stdPhase=Phase;minPhase=Phase;
    %6 points in each freq, jump from the highest one, 5 jumps atmost
    stepStd=zeros(10,6);
    for m=1:10
        stepStd(m,1)=std(Phase(m,:));
        for jump=1:5
            phiTop=find(Phase(m,:)==max(Phase(m,:)));
            stdPhase(m,phiTop)=Phase(m,phiTop)-360;
            stepStd(m,jump+1)=std(stdPhase(m,:));
        end
    end
    
    for i=1:10
        [~,minStd(i)]=min(stepStd(i,:)); %get order of smallest std (best situation)
    end
    
    for n=1:10
        if minStd(n)~=1
            for rejump=1:minStd(n)-1 %jump back from the lowest point to the best situation
                phiTop=find(stdPhase(n,:)==max(stdPhase(n,:)));
                minPhase(n,phiTop)=Phase(n,phiTop)-360;
            end
        end
    end
    for i=1:10
        %         cirVar(i,1)=circ_var(deg2rad(minPhase(i,:))')*180/pi;
        cirVar(i,1)=mean(deg2rad(minPhase(i,:))')*180/pi;
    end
    
    % find n
    nstep=cell(1,10);
    for i=1:10
        nstep{i}= ceil(0.08*freSeq(i)-1):1:floor(0.22*freSeq(i));
    end
    
    A=ones(60,2);b=[];W=[];
    for i=1:10
        A(((i-1)*6+1):i*6,2)= freSeq(i);
        W=[W;SNR(i,:)'];
    end
    count=0;
    N=cell(1,10);
    for i=1:10
        N{i}=min(nstep{i}):1:max(nstep{i});
    end
    % find linear regression (least square) that fits best
    for n10=min(nstep{10}):1:max(nstep{10})
        for n9=min(nstep{9}):1:max(nstep{9})
            for n8=min(nstep{8}):1:max(nstep{8})
                for n7=min(nstep{7}):1:max(nstep{7})
                    for n6=min(nstep{6}):1:max(nstep{6})
                        for n5=min(nstep{5}):1:max(nstep{5})
                            for n4=min(nstep{4}):1:max(nstep{4})
                                for n3=min(nstep{3}):1:max(nstep{3})
                                    for n2=min(nstep{2}):1:max(nstep{2})
                                        for n1=min(nstep{1}):1:max(nstep{1})
                                            N360=[  ...
                                                n1 n1 n1 n1 n1 n1; ...
                                                n2 n2 n2 n2 n2 n2;...
                                                n3 n3 n3 n3 n3 n3;...
                                                n4 n4 n4 n4 n4 n4;...
                                                n5 n5 n5 n5 n5 n5;...
                                                n6 n6 n6 n6 n6 n6;...
                                                n7 n7 n7 n7 n7 n7;...
                                                n8 n8 n8 n8 n8 n8;...
                                                n9 n9 n9 n9 n9 n9;...
                                                n10 n10 n10 n10 n10 n10];
                                            b= minPhase-360*N360;
                                            B= reshape(b',[],1);
                                            [Xp,~,MSE]=lscov(A,B,W);
                                            slope_tmp=Xp(2);
                                            latency_tmp=(slope_tmp/-360)*1000;
                                            if Xp(2)<-28.8 && Xp(2)>-79.2 && latency_tmp>79 && latency_tmp <219 % slope due to latency 80ms-220ms
                                                count=count+1;
                                                yaxis(count,:)=B;
                                                LR_Con(count,:)=[n1 n2 n3 n4 n5 n6 n7 n8 n9 n10];% record step moves in each linear regression conditions
                                                slope(count,1)=Xp(2);
                                                intercept(count,1)=Xp(1);
                                                Err(count,1)=MSE;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    Err=round(Err,6);
    num=find(Err==min(Err));
    sizeN=size(num);
    k=sizeN(1);
    PhaReal=[];PhaMean=[];
    figure('NumberTitle','on')
    for i=1:k
        n_fit(i,:)=LR_Con(num(i),:);
        for m=1:10
            PhaReal_tmp(m,:)=minPhase(m,:)-n_fit(i,m)*360;
            PhaMean_tmp(m,1)= mean(PhaReal_tmp(m,:),2);
        end
        PhaReal=[PhaReal;PhaReal_tmp];
        PhaMean=[PhaMean;PhaMean_tmp];
        Xp_best(i)=slope(num(i));
        intercept_best(i)=intercept(num(i));
        LATENCY(i)=(Xp_best(i)/-360)*1000;
        scatter(A(:,2),yaxis(num(i),:)','bo');
        hold on;
        plot(freSeq, PhaMean_tmp,'r*');
        str=[num2str(PhaMean_tmp)];
        text(freSeq, PhaMean_tmp,cellstr(str));
        hold on;
        grid on;
        fitLine(i,:)=freSeq*Xp_best(i)+intercept_best(i);
        plot(freSeq,fitLine(i,:),'g-');
    end
    hold off;
    pause(2)
    latency=(Xp_best(1)/-360)*1000;
    title(['slope:',num2str( Xp_best(1)),' latency',num2str(latency)]);
    xlabel('Frequency (Hz)')
    ylabel('phase(rad)')
    all_sub_latency(fold_num)=latency;
    toc
end