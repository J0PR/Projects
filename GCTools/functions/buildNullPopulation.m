function nullPop=buildNullPopulation(data,numSurrogates,nullPopMethod)
if numSurrogates==0
    nullPop=[];
    return;
end
if strcmp(nullPopMethod,'AAFT')
    nullPop=nullPopulationAFFT(data,numSurrogates);
elseif strcmp(nullPopMethod,'phaseran')
    nullPop=nullPopulationPhase(data,numSurrogates);
elseif strcmp(nullPopMethod,'IMFsurrogates')
    nullPop=IMFsurrogates(data,numSurrogates);
elseif strcmp(nullPopMethod,'TV-IMFsurrogates')
    
    opts.mix_phase_intervals=1;
    opts.mix_phase_patches=0;
    
    %no amp mixing for non-stat null hypotesis
    opts.mix_amp_patches=0;
    opts.mix_amp_intervals=0;
    opts.mix_amp_whole=0;
    opts.smooth_amp=0;
    
    
    nullPop=IMFsurrogates(data,numSurrogates,opts);
end
