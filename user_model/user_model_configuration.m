function values=user_model_configuration()
%USER_MODEL_CONFIGURATION GR4JB example's source defaults and UI metadata.
% When replacing the example, update this catalog together with both backends.
    symbols={'T_sm';'eps_m';'eps_s';'kappa';'bg';'bR';'mts';'eta';'tau';'rho'};
    names={'Snow/rain smoothing width';'Smooth minimum/maximum width';'Storage regularization width';'Percolation shape constant';'Groundwater exchange exponent';'Routing outflow exponent';'Model time-step factor';'Unit-hydrograph shape exponent';'Unused routing coefficient';'Unused storage-transition smoothing'};
    units={[char(176) 'C'];'mixed units';'mm';'-';'-';'-';'sample';'-';'-';'-'};
    defaults=[1;1e-6;1e-12;4/9;3.5;5;1;2.5;2.5;.01];
    editable=[true(6,1);false;true;false;false];lower=[0;0;0;-eps;0;1;0;-1;-Inf;0];
    references={'';'';'';'Perrin et al. (2003)';'Mathevet (2005)';'Perrin et al. (2003)';'';'';'';''};
    values=table(names,symbols,units,defaults,editable,lower,references,'VariableNames',{'Name','Symbol','Unit','Default','Editable','Lower','Reference'});
end
