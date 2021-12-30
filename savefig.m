function savefig(h,name,path)
% Attention: there should be no nullspace in either 'path' or 'home'!!!

home = pwd;
% if path is specified, the figures are saved at this directory,
% otherwise they are saved in the current directory
if nargin==3, eval(['cd ' path]); end

% save figures
saveas(h, name, 'fig');
saveas(h, name, 'jpg');
print(h, '-depsc2', name);

% go back to the current directory
if nargin==3, eval(['cd ' home]); end