function[NaNstatus] = hasNaN(Data)
%% Checks if an N dimensional dataset contains any NaNs. Returns a boolean true if any NaNs are found.

NaNstatus = any( isnan(Data(:)) );
