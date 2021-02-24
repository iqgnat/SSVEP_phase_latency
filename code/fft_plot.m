function [ frequency,fft_result ] = fft_plot( data,Fs,varargin )
% Calculate or plot directly fft results of data.
%
% [ frequency,fft_result ] = fft_plot( data,Fs,'plot' )
%
% inputs:
%   (1) data: data used to analysis. one row -> one data
%   (2) Fs: sample frequency
%   (3) 'plot': veriable input. if there is not this input, fft results will not be
%   ploted
% output:
%   (1) freqeuncy: frequency corresponding to the fft results
%   (2) fft_result: fft results

if nargin<2
    error('data and Fs must be given');
elseif nargin==2
    for k=1:size(data,1)
        size_data=size(data(k,:));
        if size_data(1)~=1 && size_data(2)~=1
            error('the length or the number of rows must be one.');
        end
        data(k,:)=detrend(data(k,:));
        L=length(data(k,:));
        NFFT=10*Fs;
%         NFFT=2^nextpow2(L);
        fft_result_temp=fft(data(k,:),NFFT)/length(data(k,:));
        fft_result(k,:)=fft_result_temp(k,1:NFFT/2+1);
        frequency(k,:)=Fs/2*linspace(0,1,NFFT/2+1);
    end
elseif nargin==3
    figure;
    title('FFT')
    for k=1:size(data,1)
        if strcmp(varargin,'plot')
             size_data=size(data(k,:));
            if size_data(1)~=1 && size_data(2)~=1
                error('the length or the number of rows must be one.');
            end
            data(k,:)=detrend(data(k,:));
            L=length(data(k,:));
            NFFT=2^nextpow2(L);
            fft_result_temp=fft(data(k,:),NFFT)/length(data(k,:));
            frequency(k,:)=Fs/2*linspace(0,1,NFFT/2+1);
            fft_result(k,:)=fft_result_temp(1:NFFT/2+1);
            subplot(size(data,1),1,k);
            plot(frequency(k,:),2*abs(fft_result(k,:)));
            xlabel('Frequency (Hz)','Fontsize',16)
            ylabel('Amplitude','Fontsize',16)
        else
            error('variable input must be ''plot''');
        end
    end
elseif nargin>=3
    error('Too much inputs')
end
end


