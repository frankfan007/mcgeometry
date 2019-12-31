function [fullmsg]=ProgressBar(Percent_Done)

% prints progress bar
% input:
%   Percent_Done: percentage of code done

CW_size=matlab.desktop.commandwindow.size;
CW_size=CW_size(1);
CW_size=CW_size-9;

N_bars=floor(Percent_Done/100*CW_size);
bars=repmat('=',[1,N_bars-1]);
blanks=repmat(' ',[1,CW_size-N_bars]);
fullmsg=sprintf('[%s>%s] %03.1f%%\n',bars,blanks,Percent_Done);
% delmsg=sprintf('>%s] %03.1f%%\n',blanks,Percent_Done);
end