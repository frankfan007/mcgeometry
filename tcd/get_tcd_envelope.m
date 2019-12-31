function [output_cell,output_answer]=get_tcd_envelope(t_nirs,tcd,varargin)

% obtains the upper and lower envelopes of the TCD signal

% input:
%   t_nirs: array with timepoints, dimension (ntimepoints,1)
%   tcd: left and right TCD signal, dimension (ntimepoints,2)

% output:
%   output_cell{1}: left TCD timepoint array
%   output_cell{2}: left TCD lower envelope array
%   output_cell{3}: left TCD upper envelope array

%   output_cell{4}: right TCD timepoint array
%   output_cell{5}: right TCD lower envelope array
%   output_cell{6}: right TCD upper envelope array

%% quality check

if isempty(varargin)
    figure(149);
    clf
    plot(t_nirs,tcd);
    ylabel('TCD signal');xlabel('seconds');
    legend('left tcd','right tcd')
    
    prompt={'Use left tcd (yes or no)?','Use right tcd (yes or no)?'};
    dlgtitle='TCD Quality Check';
    dims=[1 35];
    definput={'yes','yes'};
    opts.WindowStyle='normal';
    answer=inputdlg(prompt,dlgtitle,dims,definput,opts);
    
    output_cell=cell(6,1);
    output_answer=zeros(1,2);
    
    if strcmp(answer{1},'yes')
        output_answer(1)=1;
    end
    
    if strcmp(answer{2},'yes')
        output_answer(2)=1;
    end
else
    output_answer=varargin{1};
end

fs=1;
start_index=50;

%% left side fft

if output_answer(1)
    
    left_tcd=squeeze(tcd(:,1));
    
    Y=fft(left_tcd(~isnan(left_tcd)));
    L=length(Y);
    
    P2=abs(Y/L); P1=P2(1:L/2+1);P1(2:end-1)=2*P1(2:end-1);
    f=fs*(0:(L/2))/L;
       
    [~,index_of_max_value]=max(P1(start_index:end));
    cut_freq=f(start_index:end);
    window_length=round(1/cut_freq(index_of_max_value))+5;
    
    idx=1;
    for I=1:window_length:(floor(length(left_tcd)/window_length)*window_length - window_length)
        all_left_min(idx)=min(left_tcd(I:I+window_length));
        all_left_max(idx)=max(left_tcd(I:I+window_length));
        idx=idx+1;
    end
    
    t_left=t_nirs(1:window_length:(floor(length(left_tcd)/window_length)*window_length) - window_length);
    
    left_min=movmean(all_left_min,5);
    left_max=movmean(all_left_max,5);
    
    output_cell{1}=t_left;
    output_cell{2}=left_min;
    output_cell{3}=left_max;
end

%% right side fft 

if output_answer(2)
    
    right_tcd=squeeze(tcd(:,2));
    
    Y=fft(right_tcd(~isnan(right_tcd)));
    L=length(Y);
    
    P2=abs(Y/L); P1=P2(1:L/2+1);P1(2:end-1)=2*P1(2:end-1);
    f=fs*(0:(L/2))/L;
        
    [~,index_of_max_value]=max(P1(start_index:end));
    cut_freq=f(start_index:end);
    window_length=round(1/cut_freq(index_of_max_value))+5;
    
    idx=1;
    for I=1:window_length:(floor(length(right_tcd)/window_length)*window_length- window_length)
        all_right_min(idx)=min(right_tcd(I:I+window_length));
        all_right_max(idx)=max(right_tcd(I:I+window_length));
        idx=idx+1;
    end
    
    t_right=t_nirs(1:window_length:(floor(length(right_tcd)/window_length)*window_length) - window_length);
    
    right_min=movmean(all_right_min,5);
    right_max=movmean(all_right_max,5);
    
    output_cell{4}=t_right;
    output_cell{5}=right_min;
    output_cell{6}=right_max;
end

%% plot result

figure(149);
clf
plot(t_nirs,tcd);
hold on

plot(output_cell{1},output_cell{2},'k-','LineWidth',1)
plot(output_cell{1},output_cell{3},'k-','LineWidth',1)

plot(output_cell{4},output_cell{5},'k-','LineWidth',1)
plot(output_cell{4},output_cell{6},'k-','LineWidth',1)

ylabel('TCD signal');xlabel('seconds');
legend('left tcd','right tcd')
hold off

