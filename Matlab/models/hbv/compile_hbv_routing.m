function file = compile_hbv_routing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPILE_HBV_ROUTING Build the thin MATLAB adapter to the shared HBV router.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isdeployed
        error('hbv_route:MissingMex', ...
            'The packaged shared HBV routing MEX is missing.');
    end
    native = sage_native_root();
    output = sage_native_build_dir();
    gate = fullfile(native,'adapters','matlab','hbv', ...
        'hbv_route_mex.cpp');
    core = fullfile(native,'models','hbv');
    mex('-R2018a','-O',['-I' core],gate, ...
        '-outdir',output,'-output','hbv_route_mex');
    file = fullfile(output,['hbv_route_mex.' mexext]);
    rehash;
end
