% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/build_cartplus_mex.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function build_cartplus_mex()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD_CARTPLUS_MEX Build the MATLAB gateway for the native CaRT+ kernel.
%
% SYNOPSIS:
% build_cartplus_mex(...)
%
% INPUT ARGUMENTS:
% None.
%
% OUTPUT ARGUMENTS:
% No return value; performs the documented build or validation.
%
% NOTES:
% BUILD_CARTPLUS_MEX Compile the native CaRT+ solver next to CARTPLUS. Run
% mex -setup C++ first to select a supported MATLAB compiler.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(fileparts(here)));
    native = fullfile(root,'Native');
    kernel = fullfile(native,'models','CaRTPlus');
    gateway = fullfile(native,'adapters','matlab','cartplus', ...
        'crr_cartplus_mex.cpp');
    mex('-R2017b',['-I' kernel],'-outdir',here,gateway);
end
