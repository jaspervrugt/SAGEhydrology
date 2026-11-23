% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cart/build_cart_mex.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function build_cart_mex(outputDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD_CART_MEX Build the MATLAB gateway for the native CaRT kernel.
%
% SYNOPSIS:
% build_cart_mex(...)
%
% INPUT ARGUMENTS:
% outputDir : Destination directory for the compiled MEX file.
%
% OUTPUT ARGUMENTS:
% No return value; performs the documented build or validation.
%
% NOTES:
% BUILD_CART_MEX Compile the daily CaRT reference update for MATLAB case 4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(fileparts(here)));
    source = fullfile(root,'Native','models','cart','crr_cart_mex.cpp');
    if nargin < 1
        outputDir = here;
    end
    compiler = mex.getCompilerConfigurations('C++','Selected');
    flags = {'-O'};
    if ~isempty(compiler) && contains(compiler.Name,'MinGW')
        flags{end+1} = 'CXXOPTIMFLAGS=-O3 -DNDEBUG';
    end
    mex(flags{:},'-R2017b','-outdir',outputDir,source);
end
