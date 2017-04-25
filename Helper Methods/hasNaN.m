function[NaNstatus] = hasNaN(Data)
%% Checks if an N dimensional dataset contains any NaNs. Returns a boolean true if any NaNs are found.

datadim = ndims(Data);

reduction = Data;
for k = 1:datadim
    reduction = any(isnan(reduction));
end

if reduction == 1
    NaNstatus = true;
elseif reduction == 0
    NaNstatus = false;
else
    error('NaNcheck not functioning correctly');
end