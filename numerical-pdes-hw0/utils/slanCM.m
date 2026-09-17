function cmap = slanCM(name)
%SLANCM Wrapper for external colormap package "slanCM".
%
% cmap =SLANCM(name) loads a colormap from the external slanCM package
% installed in the user's MATLAB environment. This wrapper is included for
% reproducibility and to provide a consistent interface for the plotting
% scripts in this repository.
%
% Notes:
%   - The actual colormap definitions are part of the external slanCM
%     package and are NOT included in this repository.
%   - Users must install the slanCM package separately.
%   - This wrapper simply forwards the request to the installed package.
%
% Example"
%   colormap(slanCM('viridis'));
%
% Input:
%   name - (char or string) name of the colormap to load
%
% Output:
%   cmap - Nx3 colormap array

try
    cmap = feval(name);
catch
    warning('slanCM:NotInstalled', ...
        ['The slanCM package is not installed. ', ...
        'Using MATLAB default colormap instead.']);
    cmap = parula(256);
end
end