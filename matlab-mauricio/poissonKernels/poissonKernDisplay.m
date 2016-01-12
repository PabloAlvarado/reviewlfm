function poissonKernDisplay(kern, spacing)

% HEATKERNDISPLAY Display parameters of the POISSON kernel.
% FORMAT
% DESC displays the parameters of the Poisson kernel and the kernel type to 
% the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : poissonKernParamInit, modelDisplay, kernDisplay
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
fprintf('Number of terms in the series %2.4f\n', kern.nTerms)
fprintf(spacing);
fprintf('Space input domain length %2.4f\n', kern.lengthX)
fprintf(spacing);
fprintf('POISSON inverse width space X: %2.4f (length scale %2.4f)\n', ...
    kern.inverseWidthTime, 1/sqrt(kern.inverseWidthTime));
fprintf(spacing);
fprintf('POISSON inverse width space Y: %2.4f (length scale %2.4f)\n', ...
    kern.inverseWidthSpace, 1/sqrt(kern.inverseWidthSpace));
fprintf(spacing);
fprintf('POISSON sensitivity: %2.4f\n', kern.sensitivity)


