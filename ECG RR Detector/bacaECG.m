clear;
clc;
close all;
addpath("D:\\GoogleDriveAlvin\\MyResearchGroup\\SupineStandingECG\\");

ECGdata = dlmread("2S.txt");

%size_data = size(data,1);
%low_threshold = 900;
%pvc_point = zeros(0,1);

%i = 1;
%while i<650000
%    if(data(i,1)<900)
%        pvc_point = [pvc_point;i];
%        i = i+100;
%    end
%    i = i+1;
%end

% % ##############################################################
% % File created by Faizal Mahananto. 
% % Use this file to detect peak over clean ECG data. yes, work only on
% clean ECG data. 
% % Input : ECG data. Single column with row describe signal in 1 ms
% sampling interval. Becarefull, this program temporarily work on 1 ms
% sampling interval
% % max and mean amplitude and magic threshold obtained manualy from data
% % ##############################################################
%function rrlist = rrDetector(ECGdata)

data = ECGdata;

% =================== find threshold ================
max_amplitude = max(data)
mean_amplitude = mean(data)
% this is the magic threshold using that defined equation
% change the value of controler parameter to raise or down the level of
% threshold. Decrease the numerator (numerator/denumerator) to lower the level.
magic_threshold_point = (max_amplitude + mean_amplitude)*8/16 %(threshold controler parameter)
% ===================================================

% please consider sample interval. these below set as 1000ms sample interval 
window_threshold = floor(309/(1000/360)); % window is applied in the time series to crop time series and search the max value in each window
x_threshold_min = floor(300/(1000/360)); % to prevent detection of RR interval below 300ms.
x_threshold_max = floor(1300/(1000/360)); % max rr interval acceptable 900-1300
size_data = size(data,1);

index_bin = zeros(size_data,1);
previous_index_R_point = 0;
for i=1:window_threshold:size_data-window_threshold
    [C,I]=max(data(i:i+window_threshold));
    if C>=magic_threshold_point
        index_R_point = (i-1)+I;
        if index_R_point-previous_index_R_point > x_threshold_min
            index_bin(index_R_point) = C;
            previous_index_R_point = index_R_point;
        end
    end
end



prev_index = 0;
rrlist = zeros(0,2);
for i=1:size_data

    if index_bin(i) > 0
        if size(rrlist,1) == 0
           rrlist = [rrlist ; [0 0]];
           prev_index = i;
        else
           next_rrlist = i-prev_index;
           if next_rrlist < x_threshold_max 
               rrlist = [rrlist ; [next_rrlist i]];
           elseif next_rrlist < 2*x_threshold_max% is used to divide one long rr interval because of missed R wave into two lists
               next_rrlist = next_rrlist/2;
               rrlist = [rrlist ; [next_rrlist i]];
               rrlist = [rrlist ; [next_rrlist i]];
           else
               next_rrlist = next_rrlist/3;% is used to divide one long rr interval because of missed R wave into three lists
               rrlist = [rrlist ; [next_rrlist i]];
               rrlist = [rrlist ; [next_rrlist i]];
               rrlist = [rrlist ; [next_rrlist i]];
           end
           prev_index = i;
        end
    end
end
mak_rrlist = max(rrlist)