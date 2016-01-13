function rbfpKernDisplay(kern, spacing)

% RBFHKERNDISPLAY Display parameters of the RBFP kernel.
% FORMAT
% DESC displays the parameters of the radial basis function Poisson
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfpKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBFH inverse width space X: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthX, 1/sqrt(kern.inverseWidthX));
fprintf(spacing);
fprintf('RBFH inverse width space Y: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthY, 1/sqrt(kern.inverseWidthY));

