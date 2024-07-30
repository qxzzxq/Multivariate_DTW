function dat_processed = PreProcessDat(dat, pre)
%         dat  : data set in form of cell to be used which contains all the
%                time series data of an attribute. 
%         pre  : pre-processing option 
%               (1: None, 2: Z norm, 3: 0/1(min-max), 4: -1/1, 5: Decimal norm, 6: Tanh norm, 7: Sigmoid norm )
if pre == 2 
    dat_processed = cellfun(@pre_zNorm, dat, 'UniformOutput',false);
elseif pre == 3
    dat_processed = cellfun(@pre_zeroOneNorm, dat, 'UniformOutput',false);
elseif pre == 4
    dat_processed = cellfun(@pre_absOneNorm, dat, 'UniformOutput',false);
elseif pre == 5
    dat_processed = cellfun(@pre_decimalNorm, dat, 'UniformOutput',false);
elseif pre == 6
    dat_processed = cellfun(@pre_tanhNorm, dat, 'UniformOutput',false);
elseif pre == 7
    dat_processed = cellfun(@pre_sigmoidNorm, dat, 'UniformOutput',false);
else
    dat_processed = dat;
end
