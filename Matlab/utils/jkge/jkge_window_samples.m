function n = jkge_window_samples(loss)
% Prepared configurations retain user days and expose converted sample
% count.
if isfield(loss,'n_win_samples'), n=loss.n_win_samples; else, n= ...
    loss.n_win; end
end
