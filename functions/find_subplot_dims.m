function [x,y]=find_subplot_dims(num_plots)

% suggests possible length and width of a subplot array in a figure

% input:
%   num_plots: number of total subplots
%
% output:
%   x: length of subplot array
%   y: width of subplot array

%% counting for 1

if num_plots==1
    x=1; y=1; return
end

%% finding number of subplots

div_numplots=divisors(num_plots);
all_comb=combnk(div_numplots,2);
temp=all_comb(:,1).*all_comb(:,2);
possible_combs=all_comb(temp==num_plots,:);
[~,I]=min(abs(possible_combs(:,2)-possible_combs(:,1)));
x=min(possible_combs(I,:));
y=max(possible_combs(I,:));

while abs(y-x)>5
    num_plots=num_plots+1;
    div_numplots=divisors(num_plots);
    all_comb=combnk(div_numplots,2);
    temp=all_comb(:,1).*all_comb(:,2);
    possible_combs=all_comb(temp==num_plots,:);
    [~,I]=min(abs(possible_combs(:,2)-possible_combs(:,1)));
    x=min(possible_combs(I,:));
    y=max(possible_combs(I,:));
end