function varargout = usvseg09r2_plus(action)
% usvseg09r2_plus
%
% by R.O. Tachibana (rtachi@gmail.com)
% May 23, 2020

% update:
%   - changed structure of global variable
%       gv.gui / gv.dat / gv.prm
%   - output duration is now in milliseconds but not seconds
%   - added "meanfreq" and "cvfreq" as output feature
%   - deleted "whitening" and "hanntaperspec" functions for comparison
%   - parallel computeing (parfor) is now working only for batch processing

% modified by J. Matsumoto (jm@med.u-toyama.ac.jp)
% May 8, 2021

% updeate:
%   - output subsegments for usvcam 
%   - mark subsegments totally within margin (to exclude noise detection)
%   - changed liftercutoff: 3 -> 6
%   - added audible wav file output
%   - added gui and related code to process usvcam data

% global variable
global gv;

% Initiation ------------------------------------------------------
if nargin==0
    % version
    usevsegver = 'USVSEG ver 0.9 (rev 2)';
    fftsize = 512;
    % load parameters
    prmfile = [fileparts(mfilename('fullpath')) filesep 'usvseg_prm.mat'];
    if exist(prmfile)==2, l=load(prmfile); gv.prm = l.prm;
    else % default values
        prm.timestep = 0.0005; prm.freqmin = 30000; prm.freqmax = 100000;
        prm.threshval = 3.0;   prm.durmin = 0.005;  prm.durmax = 1.0;
        prm.gapmin = 0.030;    prm.margin = 0.015;  prm.wavfileoutput = 0;
        prm.imageoutput = 1;   prm.imagetype = 1;   prm.traceoutput = 1;
        prm.readsize = 2;     prm.mapL = 1000;     prm.mapH = 6000;
        prm.overwriteresult = 0;
        gv.prm = prm;
    end
    gv.usevsegver = usevsegver;
    gv.prm.fftsize = fftsize; % FFT size = window size
    gv.dat.pth = cd;
    % make GUI
    makegui;
    cla(gv.gui.haxes(1)); cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    return
end
% open ----------------------------------------------------------------
if strcmp(action,'open')
    % clear axes
    cla(gv.gui.haxes(1)); cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    % unable
    set([gv.gui.hpush_thrs gv.gui.hpush_dtct gv.gui.hpush_save gv.gui.hpush_sgmt gv.gui.hpush_play gv.gui.hpush_swav],'Enable','Off');
    % fetch parameters from GUI
    fetchparams;
    timestep = gv.prm.timestep;
    readsize = gv.prm.readsize;
    fftsize = gv.prm.fftsize;
    % get file
    [fn,pth] = uigetfile([gv.dat.pth filesep '*.wav']);
    % show file name
    set(gv.gui.htext_fnam,'string',['File: ' pth fn]);
    busytoggle(1);
    % read wav
    ai = audioinfo([pth fn]);
    fs = ai.SampleRate;
    wav = audioread([pth fn],[1 min(readsize*fs,ai.TotalSamples)]);
    % show fs
    set(gv.gui.htext_fsvl,'string',num2str(fs/1000));
    % multitaper spectrogram
    step = round(timestep*fs);
    mtsp = multitaperspec(wav,fs,fftsize,timestep,0);
    mtsimrng = median(mtsp(:)) + [0 40];
    fvec = [0; (1:(fftsize/2))'/fftsize*fs];
    tvec = ((0:(size(mtsp,2)-1))*step+fftsize/2)/fs;
    % default view range : 2 seconds
    xl = [0 2];
    yl = [0 max(fvec)/1000];
    % save in global variable
    gv.dat.fn = fn; gv.dat.pth = pth;
    gv.dat.fs = fs; gv.dat.wav = wav; gv.dat.step = step; 
    gv.dat.mtsp = mtsp; gv.dat.mtsimrng = mtsimrng;
    gv.dat.fvec = fvec; gv.dat.tvec = tvec;
    gv.dat.xl = xl; gv.dat.yl = yl;
    % draw mtsp
    drawmtsp;
    % enable
    set(gv.gui.hpush_flat,'Enable','On');
    busytoggle(0);
    return;
end
% flattening ------------------------------------------------------------
if strcmp(action,'flat')
    busytoggle(1);
    % clear axes
    cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    % flattening
    fltnd = flattening(gv.dat.mtsp);
    fltndimrng = [0 30];
    % save
    gv.dat.fltnd = fltnd;
    gv.dat.fltndimrng = fltndimrng;
    % draw fltnd
    drawfltnd;
    % enable
    set(gv.gui.hpush_thrs,'Enable','On');
    busytoggle(0);
    return;
end
% thresholding ------------------------------------------------------------
if strcmp(action,'thrs')
    % clear axes
    cla(gv.gui.haxes(3));
    % fetch
    fetchparams;
    freqmin = gv.prm.freqmin; freqmax = gv.prm.freqmax; threshval = gv.prm.threshval; 
    fs = gv.dat.fs; fvec = gv.dat.fvec; tvec = gv.dat.tvec; fltnd = gv.dat.fltnd;
    busytoggle(1);    
    % threshold calculation with threshval*sigma (SD) of background noise 
    thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval);
    % thresholding
    thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh);
    % draw
    set(0,'CurrentFigure',gv.gui.hfig);
    set(gv.gui.hfig,'CurrentAxes',gv.gui.haxes(3));
    imagesc(tvec,fvec/1000,thrshd); axis xy
    colormap(gca,gray);
    set(gca,'tickdir','out'); box off;
    ylabel('Frequency (kHz)');
    set(gv.gui.haxes,'xlim',gv.dat.xl,'ylim',gv.dat.yl);
    drawnow;
    % enable
    set(gv.gui.hpush_dtct,'Enable','On');
    % save
    gv.dat.thrshd = thrshd; gv.dat.thresh = thresh;
    busytoggle(0);
    return;
end
% detection ------------------------------------------------------------
if strcmp(action,'dtct')
    busytoggle(1);
    % clear axes
    cla(gv.gui.haxes(3));   
    % fetch
    fetchparams;
    fs = gv.dat.fs; timestep = gv.prm.timestep;
    freqmin = gv.prm.freqmin; freqmax = gv.prm.freqmax;
    durmin = gv.prm.durmin; durmax = gv.prm.durmax;
    gapmin = gv.prm.gapmin; margin = gv.prm.margin;
    mtsp = gv.dat.mtsp; fltnd = gv.dat.fltnd; thrshd = gv.dat.thrshd;
    % onset/offset detection
    [onoffset,onoffsetm] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
    % peak tracking
    [freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin);
    % save
    gv.dat.onoffset = onoffset; gv.dat.onoffsetm = onoffsetm;
    gv.dat.freqtrace = freqtrace; gv.dat.amptrace = amptrace;
    gv.dat.maxampval = maxampval; gv.dat.maxfreq = maxfreq; gv.dat.maxampidx = maxampidx;
    gv.dat.meanfreq = meanfreq; gv.dat.cvfreq = cvfreq;
    % draw trace
    drawtrace;
    % enable
    set(gv.gui.hpush_save,'Enable','On');
    set(gv.gui.hpush_sgmt,'Enable','On');
    set(gv.gui.hpush_play,'Enable','On');
    busytoggle(0);
    return;
end
% save ----------------------------------------------------------------
if strcmp(action,'save')
    busytoggle(0);
    % fetch
    onoffset = gv.dat.onoffset;
    maxampval = gv.dat.maxampval; maxfreq = gv.dat.maxfreq;
    meanfreq = gv.dat.meanfreq; cvfreq = gv.dat.cvfreq;
    dur = diff(onoffset,[],2);
    [~,prefix,~] = fileparts(gv.dat.fn);
    % save CSV
    [fn,pth] = uiputfile([gv.dat.pth filesep prefix '_dat.csv'],'file name');
    fid = fopen([pth fn],'wt');
    fprintf(fid,'#,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq\n');
    if ~isempty(onoffset)
        for n=1:size(onoffset,1)
            fprintf(fid,'%d,%.04f,%.04f,%.1f,%.03f,%.02f,%.03f,%.04f\n',n,onoffset(n,1),onoffset(n,2),dur(n)*1000,maxfreq(n)/1000,maxampval(n),meanfreq(n)/1000,cvfreq(n));
        end
    end
    fclose(fid);
    busytoggle(1);
    return;
end
% segment -------------------------------------------------------------
if strcmp(action,'sgmt')
    busytoggle(0);
    % fetch
    fetchparams;
    margin = gv.prm.margin; timestep = gv.prm.timestep;
    wavflg = gv.prm.wavfileoutput; imgflg = gv.prm.imageoutput; 
    fltflg = gv.prm.imagetype; trcflg = gv.prm.traceoutput;
    wav = gv.dat.wav; fs = gv.dat.fs; mtsp = gv.dat.mtsp; fltnd = gv.dat.fltnd;
    tvec = gv.dat.tvec; onoffset = gv.dat.onoffset;
    freqtrace = gv.dat.freqtrace;amptrace = gv.dat.amptrace;
    fltndimrng = gv.dat.fltndimrng; mtsimrng = gv.dat.mtsimrng;
    % file name
    [~,prefix,~] = fileparts(gv.dat.fn);
    outp = uigetdir(gv.dat.pth);
    % output
    startid = 1;
    if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
    else,         inputimg = mtsp;   imrng = mtsimrng; end
    segfun(startid,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,tvec,amptrace,freqtrace,wavflg,imgflg,trcflg);
    busytoggle(0);
    return;
end
% play  ----------------------------------------------------------------
if strcmp(action,'play')
    playfs = 44100;
    % fetch
    fetchparams;
    mapL = gv.prm.mapL; mapH = gv.prm.mapH;
    fs = gv.dat.fs; tvec = gv.dat.tvec;
    freqtrace = gv.dat.freqtrace; amptrace = gv.dat.amptrace;
    % audible sound synthesis
    synthsnd = soundsynthesis(freqtrace(:,1),amptrace(:,1),tvec,fs,playfs,[mapL mapH]);
    % get range
    rng = max(1,round(gv.dat.xl(1)*playfs)):min(length(synthsnd),round(gv.dat.xl(2)*playfs));
    % play
    ap = audioplayer(synthsnd(rng),playfs);
    playblocking(ap);
    delete(ap);
    set(gv.gui.hpush_swav,'Enable','On');
    % save
    gv.dat.playfs = playfs;
    gv.dat.synthsnd = synthsnd;
    return;
end
% save wav  ----------------------------------------------------------------
if strcmp(action,'swav')
    newfn = [gv.dat.fn(1:end-4) '_syn.wav'];
    [f,p] = uiputfile([gv.dat.pth newfn]);
    audiowrite([p f],gv.dat.synthsnd,gv.dat.playfs);
    return;
end
% long ----------------------------------------------------------------
if strcmp(action,'long')
    % fetch parameter
    fetchparams;
    timestep = gv.prm.timestep; margin =  gv.prm.margin; durmax = gv.prm.durmax;
    wavflg = gv.prm.wavfileoutput; imgflg = gv.prm.imageoutput;
    fltflg = gv.prm.imagetype; trcflg = gv.prm.traceoutput;
    readsize = gv.prm.readsize; fftsize = gv.prm.fftsize;
    % clear axes
    cla(gv.gui.haxes(1)); cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    % disable
    set(gv.gui.hpush_open,'enable','on');
    set(gv.gui.hpush_dtct,'enable','off');
    set(gv.gui.hpush_save,'enable','off');
    set(gv.gui.hpush_sgmt,'enable','off');
    set(gv.gui.hpush_long,'enable','on');
    % get long wavfile
    busytoggle(1);
    [fn,pth] = uigetfile([gv.dat.pth filesep '*.wav']);
    fp = [pth fn];
    ainfo = audioinfo(fp);
    wavsize = ainfo.TotalSamples;
    fs = ainfo.SampleRate;
    nreadsize= round(readsize*fs);
    nread = ceil(wavsize/nreadsize);
    fvec = [0; (1:(fftsize/2))'/fftsize*fs];
    step = round(timestep*fs);
    yl = [0 max(fvec)/1000];
    gv.dat.fs = fs;
    gv.dat.yl = yl;
    % show fs
    set(gv.gui.htext_fsvl,'string',num2str(fs/1000));
    % CSV setting
    [~,prefix,~] = fileparts(fn);
    [sfn,spth] = uiputfile([pth filesep prefix '_dat.csv'],'save file name');
    savefp = [spth sfn];
    fid = fopen(savefp,'wt');
    fprintf(fid,'#,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq\n');
    fclose(fid);
    % segmentation setting
    outp = uigetdir(spth,'output directory');
    % start
    prevn = 0;
    prevlast = 0;
    med = [];
    thresh = [];
    for r=1:nread
        set(gv.gui.htext_fnam,'string',['File: ' fp ' ...  (' num2str(r) '/' num2str(nread) ' blocks)']);
        busytoggle(1);
        % read
        rng = [prevlast+1 min(r*nreadsize,wavsize)];
        if diff(rng)<fftsize*2, break; end
        [wav,fs] = audioread(fp,rng);
        % process
        [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,gv.prm,med,thresh);
        busytoggle(0);
        dur = diff(onoffset,[],2);
        ronoffset = onoffset+(prevlast+1)/fs;
        ronoffsetm = onoffsetm+(prevlast+1)/fs;
        nstep = size(mtsp,2);
        tvec = ((0:(nstep-1))*step+fftsize/2)/fs;
        rtvec = tvec+(prevlast+1)/fs;
        mtsimrng = median(mtsp(:)) + [0 40];
        fltndimrng = [0 30];
        % save to global 
        gv.dat.tvec = rtvec; gv.dat.fvec = fvec; gv.dat.mtsp = mtsp; gv.dat.fltnd = fltnd;
        gv.dat.onoffset = ronoffset; gv.dat.onoffsetm = ronoffsetm;
        gv.dat.freqtrace = freqtrace; gv.dat.amptrace = amptrace;
        gv.dat.maxampval = maxampval; gv.dat.maxampidx = maxampidx; gv.dat.maxfreq = maxfreq;
        gv.dat.meanfreq = meanfreq; gv.dat.cvfreq = cvfreq;
        gv.dat.mtsimrng = mtsimrng; gv.dat.fltndimrng = fltndimrng;
        gv.dat.xl = [min(rtvec) max(rtvec)];
        % draw MTSP
        drawmtsp;
        % draw Flattened
        drawfltnd;
        % draw trace
        drawtrace;
        % save CSV
        fid = fopen(savefp,'at');
        if ~isempty(ronoffset)
            for n=1:size(ronoffset,1)
                fprintf(fid,'%d,%.04f,%.04f,%.1f,%.03f,%.02f,%.03f,%.04f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n)*1000,maxfreq(n)/1000,maxampval(n),meanfreq(n)/1000,cvfreq(n));
            end
        end
        fclose(fid);
        % segmentation
        if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
        else,         inputimg = mtsp;   imrng = mtsimrng; end
        segfun(prevn+1,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,rtvec,amptrace,freqtrace,wavflg,imgflg,trcflg)
        % calc
        prevn = size(onoffset,1)+prevn;
        if contflg==1
            if isempty(onoffset)
                prevlast = rng(2) - durmax * fs;
            else
                prevlast = round(onoffset(end,2)*fs)+prevlast;
            end
        else
            prevlast = rng(2);
        end        
    end
    
    %%%%%% output audible sound file; added by JM May,7,2021
    playfs = 44100;
    mapL = gv.prm.mapL; mapH = gv.prm.mapH;
    L = dir([outp filesep '*.ori.csv']);
    data = [];
    for i=1:length(L)
       d = dlmread([outp filesep L(i).name]);
       data = [data; d(:,1:3)];
    end
    if isempty(data)
        synthsnd = zeros(playfs,1);
    else
        data = sortrows(data, 1);
        [~,I,~] = unique(data(:,1));
        data = data(I,:);
        synthsnd = soundsynthesis(data(:,3),data(:,2),data(:,1)',fs,playfs,[mapL mapH]);
    end
    [p, n, ~] = fileparts(fp);
    fp_synth = [p filesep n '.audible.wav'];
    audiowrite(fp_synth, synthsnd, playfs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    msgbox('Done!');
    return;
end
% folder (multi/single usvcam files) ----added by JM May,8 2021------------------------------------------
if strcmp(action,'usvcam')
    % get reading path and files
    pth = uigetdir(gv.dat.pth);
    gv.dat.pth = pth;
    
    L = dir([pth filesep '*.wav']);
    
    flist = {};
    
    if isempty(L)
        LL = dir([pth filesep '*']);
        for i=3:length(LL)
            pth2 = [pth filesep LL(i).name];
            if exist(pth2, 'dir')
                flist{end+1} = pth2;
            end
        end
    else
       flist = {pth};
    end
    
    set(gv.gui.htext_fnam,'string','');
    % clear axes
    cla(gv.gui.haxes(1)); cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    % disable
    set(gv.gui.hpush_open,'enable','on');
    set(gv.gui.hpush_dtct,'enable','off');
    set(gv.gui.hpush_save,'enable','off');
    set(gv.gui.hpush_sgmt,'enable','off');
    set(gv.gui.hpush_long,'enable','on');
    % fetch parameter
    fetchparams;
    timestep = gv.prm.timestep; margin =  gv.prm.margin; durmax = gv.prm.durmax;
    wavflg = gv.prm.wavfileoutput; imgflg = gv.prm.imageoutput;
    fltflg = gv.prm.imagetype; trcflg = gv.prm.traceoutput;
    readsize = gv.prm.readsize; fftsize = gv.prm.fftsize;
    
    for i_file=1:length(flist)
        
        if ~isempty(dir([flist{i_file} filesep '*.usvseg_dat.csv']))  % results already exist
            if gv.prm.overwriteresult
                %clean previous result
                delete([flist{i_file} filesep '*.audible.wav']);
                delete([flist{i_file} filesep '*.usvseg_dat.csv']);
                delete([flist{i_file} filesep 'seg/*.csv']);
                delete([flist{i_file} filesep 'seg/*.wav']);
                delete([flist{i_file} filesep 'seg/*.jpg']);
                if exist([flist{i_file} filesep 'seg'], 'dir')
                    rmdir([flist{i_file} filesep 'seg']);
                end
            else
                continue;
            end
        end
        
        fp = dir([flist{i_file} filesep '*.wav']);
        if isempty(fp)
            continue;
        end
        
        fp = [flist{i_file} filesep fp(1).name];
        
        % get wavfile
        busytoggle(1);
        ainfo = audioinfo(fp);
        wavsize = ainfo.TotalSamples;
        fs = ainfo.SampleRate;
        nreadsize= round(readsize*fs);
        nread = ceil(wavsize/nreadsize);
        fvec = [0; (1:(fftsize/2))'/fftsize*fs];
        step = round(timestep*fs);
        yl = [0 max(fvec)/1000];
        gv.dat.fs = fs;
        gv.dat.yl = yl;
        % show fs
        set(gv.gui.htext_fsvl,'string',num2str(fs/1000));
        % CSV setting
        [~,prefix,~] = fileparts(fp);
        savefp = [flist{i_file} filesep prefix '.usvseg_dat.csv'];
        fid = fopen(savefp,'wt');
        fprintf(fid,'#,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq\n');
        fclose(fid);
        % segmentation setting
        outp = [flist{i_file} filesep 'seg'];
        mkdir(outp);
        % start
        prevn = 0;
        prevlast = 0;
        med = [];
        thresh = [];
        for r=1:nread
            set(gv.gui.htext_fnam,'string',['File: ' fp ' ...  (' num2str(r) '/' num2str(nread) ' blocks)']);
            busytoggle(1);
            % read
            rng = [prevlast+1 min(r*nreadsize,wavsize)];
            if diff(rng)<fftsize*2, break; end
            [wav,fs] = audioread(fp,rng);
            % process
            [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,gv.prm,med,thresh);
            busytoggle(0);
            dur = diff(onoffset,[],2);
            ronoffset = onoffset+(prevlast+1)/fs;
            ronoffsetm = onoffsetm+(prevlast+1)/fs;
            nstep = size(mtsp,2);
            tvec = ((0:(nstep-1))*step+fftsize/2)/fs;
            rtvec = tvec+(prevlast+1)/fs;
            mtsimrng = median(mtsp(:)) + [0 40];
            fltndimrng = [0 30];
            % save to global 
            gv.dat.tvec = rtvec; gv.dat.fvec = fvec; gv.dat.mtsp = mtsp; gv.dat.fltnd = fltnd;
            gv.dat.onoffset = ronoffset; gv.dat.onoffsetm = ronoffsetm;
            gv.dat.freqtrace = freqtrace; gv.dat.amptrace = amptrace;
            gv.dat.maxampval = maxampval; gv.dat.maxampidx = maxampidx; gv.dat.maxfreq = maxfreq;
            gv.dat.meanfreq = meanfreq; gv.dat.cvfreq = cvfreq;
            gv.dat.mtsimrng = mtsimrng; gv.dat.fltndimrng = fltndimrng;
            gv.dat.xl = [min(rtvec) max(rtvec)];
            % draw MTSP
            drawmtsp;
            % draw Flattened
            drawfltnd;
            % draw trace
            drawtrace;
            % save CSV
            fid = fopen(savefp,'at');
            if ~isempty(ronoffset)
                for n=1:size(ronoffset,1)
                    fprintf(fid,'%d,%.04f,%.04f,%.1f,%.03f,%.02f,%.03f,%.04f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n)*1000,maxfreq(n)/1000,maxampval(n),meanfreq(n)/1000,cvfreq(n));
                end
            end
            fclose(fid);
            % segmentation
            if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
            else,         inputimg = mtsp;   imrng = mtsimrng; end
            segfun(prevn+1,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,rtvec,amptrace,freqtrace,wavflg,imgflg,trcflg)
            % calc
            prevn = size(onoffset,1)+prevn;
            if contflg==1
                if isempty(onoffset)
                    prevlast = rng(2) - durmax * fs;
                else
                    prevlast = round(onoffset(end,2)*fs)+prevlast;
                end
            else
                prevlast = rng(2);
            end        
        end

        %%%%%% output audible sound file; added by JM May,7,2021
        playfs = 44100;
        mapL = gv.prm.mapL; mapH = gv.prm.mapH;
        L = dir([outp filesep '*.ori.csv']);
        data = [];
        for i=1:length(L)
           d = dlmread([outp filesep L(i).name]);
           data = [data; d(:,1:3)];
        end
        if isempty(data)
            synthsnd = zeros(playfs,1);
        else
            data = sortrows(data, 1);
            [~,I,~] = unique(data(:,1));
            data = data(I,:);
            synthsnd = soundsynthesis(data(:,3),data(:,2),data(:,1)',fs,playfs,[mapL mapH]);
        end
        [p, n, ~] = fileparts(fp);
        fp_synth = [p filesep n '.audible.wav'];
        audiowrite(fp_synth, synthsnd, playfs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    msgbox('Done!');
    return;
end
% folder (multiple files) ----------------------------------------------------------------
if strcmp(action,'fldr')
    set(gv.gui.htext_fnam,'string','');
    % clear axes
    cla(gv.gui.haxes(1)); cla(gv.gui.haxes(2)); cla(gv.gui.haxes(3));
    % disable
    set(gv.gui.hpush_open,'enable','on');
    set(gv.gui.hpush_dtct,'enable','off');
    set(gv.gui.hpush_save,'enable','off');
    set(gv.gui.hpush_sgmt,'enable','off');
    set(gv.gui.hpush_long,'enable','on');
    % fetch parameters
    fetchparams;
    timestep = gv.prm.timestep; margin =  gv.prm.margin; durmax = gv.prm.durmax;
    wavflg = gv.prm.wavfileoutput; imgflg = gv.prm.imageoutput;
    fltflg = gv.prm.imagetype; trcflg = gv.prm.traceoutput;
    readsize = gv.prm.readsize; fftsize = gv.prm.fftsize;
    % get reading path and files
    pth = uigetdir(gv.dat.pth);
    gv.dat.pth = pth;
    d = dir([pth filesep '*.wav']);
    fns = {d.name};
    nfile = length(fns);
    for f=1:nfile
        % get wavfile
        fp = [pth filesep fns{f}];
        ainfo = audioinfo(fp);
        wavsize = ainfo.TotalSamples;
        fs = ainfo.SampleRate;
        nreadsize= round(readsize*fs);
        nread = ceil(wavsize/nreadsize);
        fvec = [0; (1:(fftsize/2))'/fftsize*fs];
        step = round(timestep*fs);
        yl = [0 max(fvec)/1000];
        gv.dat.fs = fs;
        gv.dat.yl = yl;
        % show fs
        set(gv.gui.htext_fsvl,'string',num2str(fs/1000));
        % CSV setting
        [~,prefix,~] = fileparts(fns{f});
        sfn = [prefix '_dat.csv'];
        savefp = [pth filesep sfn];
        fid = fopen(savefp,'wt');
        fprintf(fid,'#,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq\n');
        fclose(fid);
        % segmentation setting
        mkdir(pth,prefix);
        outp = [pth filesep prefix];
        % start
        prevn = 0;
        prevlast = 0;
        med = [];
        thresh = [];
        for r=1:nread
            % read
            set(gv.gui.htext_fnam,'string',sprintf('File: %s  ... (%d/%d blocks) [%d/%d files]',fp,r,nread,f,nfile));
            busytoggle(1);
            rng = [prevlast+1 min(r*nreadsize,wavsize)];
            if diff(rng)<fftsize*2, break; end
            [wav,fs] = audioread(fp,rng);
            % process
            [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,gv.prm,med,thresh);
            dur = diff(onoffset,[],2);
            ronoffset = onoffset+(prevlast+1)/fs;
            ronoffsetm = onoffsetm+(prevlast+1)/fs;
            nstep = size(mtsp,2);
            tvec = ((0:(nstep-1))*step+fftsize/2)/fs;
            rtvec = tvec+(prevlast+1)/fs;
            mtsimrng = median(mtsp(:)) + [0 40];
            fltndimrng = [0 30];
            busytoggle(0);
            % save to global 
            gv.dat.tvec = rtvec; gv.dat.fvec = fvec; gv.dat.mtsp = mtsp; gv.dat.fltnd = fltnd;
            gv.dat.onoffset = ronoffset; gv.dat.onoffsetm = ronoffsetm;
            gv.dat.freqtrace = freqtrace; gv.dat.amptrace = amptrace;
            gv.dat.maxampval = maxampval; gv.dat.maxampidx = maxampidx; gv.dat.maxfreq = maxfreq;
            gv.dat.meanfreq = meanfreq; gv.dat.cvfreq = cvfreq;
            gv.dat.mtsimrng = mtsimrng; gv.dat.fltndimrng = fltndimrng;
            gv.dat.xl = [min(rtvec) max(rtvec)];
            % draw MTSP
            drawmtsp;
            % draw Flattened
            drawfltnd;
            % draw trace
            drawtrace;
            % save CSV
            fid = fopen(savefp,'at');
            if ~isempty(ronoffset)
                for n=1:size(ronoffset,1)
                    fprintf(fid,'%d,%.04f,%.04f,%.1f,%.03f,%.02f,%.03f,%.04f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n)*1000,maxfreq(n)/1000,maxampval(n),meanfreq(n)/1000,cvfreq(n));
                end
            end
            fclose(fid);
            % segmentation
            if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
            else,         inputimg = mtsp;   imrng = mtsimrng; end
            segfun(prevn+1,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,rtvec,amptrace,freqtrace,wavflg,imgflg,trcflg);
            % calc
            prevn = size(onoffset,1)+prevn;
            if contflg==1
                if isempty(onoffset)
                    prevlast = rng(2) - durmax * fs;
                else
                    prevlast = round(onoffset(end,2)*fs)+prevlast;
                end
            else
                prevlast = rng(2);
            end        
        end
    end
    msgbox('Done!');
    return;
end
% zoom  -----------------------------------------------------
if strcmp(action,'z')
    set(gcf,'Pointer','crosshair')
    waitforbuttonpress;
    p1 = get(gca,'currentpoint');
    rbbox;
    p2 = get(gca,'currentpoint');
    newxl = sort([p1(1,1) p2(1,1)]);
    newyl = sort([p1(1,2) p2(1,2)]);
    set(gcf,'Pointer','arrow')
    gv.dat.xl = newxl;
    gv.dat.yl = newyl;
    set([gv.gui.haxes],'xlim',gv.dat.xl,'ylim',gv.dat.yl);
    return;
end
% zoom in -----------------------------------------------------
if strcmp(action,'i')
    xl = get(gca,'XLim');
    center = xl(1)+diff(xl)/2;
    newxl = center + [-diff(xl)/2 diff(xl)/2]*0.75;
    gv.dat.xl = newxl;
    set([gv.gui.haxes],'xlim',gv.dat.xl);
    return;
end
% zoom out -----------------------------------------------------
if strcmp(action,'o')
    xl = get(gca,'XLim');
    center = xl(1)+diff(xl)/2;
    newxl = center + [-diff(xl)/2 diff(xl)/2]*1.25;
    gv.dat.xl = newxl;
    set([gv.gui.haxes],'xlim',gv.dat.xl);
    return;
end
% zoom in veritcal --------------------------------------------------
if strcmp(action,'u')
    yl = get(gca,'ylim');
    center = yl(1)+diff(yl)/2;
    newyl = center + [-diff(yl)/2 diff(yl)/2]*0.75;
    gv.dat.yl = newyl;
    set([gv.gui.haxes],'ylim',gv.dat.yl);
    return;
end
% zoom out vertical ---------------------------------------------------
if strcmp(action,'j')
    yl = get(gca,'ylim');
    center = yl(1)+diff(yl)/2;
    newyl = center + [-diff(yl)/2 diff(yl)/2]*1.25;
    gv.dat.yl = newyl;
    set([gv.gui.haxes],'ylim',gv.dat.yl);
    return;
end
% slide forward ----------------------------------------------
if double(action)==29 % ->
    xl = get(gca,'xlim');
    newxl = xl+diff(xl)*0.05;
    gv.dat.xl = newxl;
    set([gv.gui.haxes],'xlim',gv.dat.xl);
    return;
end
% slide back --------------------------------------------------
if double(action)==28 % <-
    xl = get(gca,'xlim');
    newxl = xl-diff(xl)*0.05;
    gv.dat.xl = newxl;
    set([gv.gui.haxes],'xlim',gv.dat.xl);
    return;
end
% slide up ----------------------------------------------
if double(action)==30 % up
    yl = get(gca,'ylim');
    newyl = yl+diff(yl)*0.1;
    gv.dat.yl = newyl;
    set([gv.gui.haxes],'ylim',gv.dat.yl);
    return;
end
% slide down --------------------------------------------------
if double(action)==31 % down
    xl = get(gca,'ylim');
    newyl = xl-diff(xl)*0.1;
    gv.dat.yl = newyl;
    set([gv.gui.haxes],'ylim',gv.dat.yl);
    return;
end
% save PDF file ------------------------------------------------
if strcmp(action,'p')
    [f,p] = uiputfile([gv.dat.pth filesep '*.pdf']);
    set(gcf,'renderer','painters','PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print([p f],'-dpdf');
    return;
end
% help --------------------------------------------------
if strcmp(action,'h')
    msg = {
    'i ... zoom in horizontal'
    'o ... zoom out horizontal'
    'u ... zoom in vertical'
    'j ... zoom out vertical'
    'right ... slide to right'
    'left ... slide to left'
    'up ... slide to upper'
    'down ... slide to lower'
    'p ... save PDF file'
    'h ... help'
    };
    msgbox(msg,'Key Command List');    
    return;
end
% for other things --------------------------------------
if length(action)==1
    return;
end
% for internal functions --------------------------------------
if length(action)>5
    varargout{1} = eval(['@' action]);
    return;
end
% ////////////////////////////////////////////////////////////////////////
function drawmtsp
global gv;

set(0,'CurrentFigure',gv.gui.hfig);
set(gv.gui.hfig,'CurrentAxes',gv.gui.haxes(1));
imagesc(gv.dat.tvec,gv.dat.fvec/1000,gv.dat.mtsp);  axis xy
set(gca,'tickdir','out'); box off
ylabel('Frequency (kHz)');
colormap(flipud(gray));
caxis(gv.dat.mtsimrng);
set(gv.gui.haxes,'xlim',gv.dat.xl,'ylim',gv.dat.yl);
drawnow;
% ////////////////////////////////////////////////////////////////////////
function drawfltnd
global gv;
% fetch
fetchparams;    
freqmin = gv.prm.freqmin;
freqmax = gv.prm.freqmax;
fvec = gv.dat.fvec;
tvec = gv.dat.tvec;
fltnd = gv.dat.fltnd;
fltndimrng = gv.dat.fltndimrng;
% draw
set(0,'CurrentFigure',gv.gui.hfig);
set(gv.gui.hfig,'CurrentAxes',gv.gui.haxes(2));
imagesc(tvec,fvec/1000,fltnd); axis xy
set(gca,'tickdir','out'); box off;
caxis(fltndimrng);
hls = line([min(tvec) max(tvec); min(tvec) max(tvec)]',[freqmax freqmax;freqmin freqmin]'/1000); set(hls,'color','r','linestyle','--');
ylabel('Frequency (kHz)');
set(gv.gui.haxes,'xlim',gv.dat.xl,'ylim',gv.dat.yl);
drawnow;
% ////////////////////////////////////////////////////////////////////////
function drawtrace
global gv;
% fetch
fvec = gv.dat.fvec;
tvec = gv.dat.tvec;
fltnd = gv.dat.fltnd;
fltndimrng = gv.dat.fltndimrng;
freqtrace = gv.dat.freqtrace;
onoffset = gv.dat.onoffset;
onoffsetm = gv.dat.onoffsetm;
maxampidx = gv.dat.maxampidx;
maxfreq = gv.dat.maxfreq;
margin = gv.prm.margin;
fs = gv.dat.fs;
% draw
set(0,'CurrentFigure',gv.gui.hfig);
set(gv.gui.hfig,'CurrentAxes',gv.gui.haxes(3));
imagesc(tvec,fvec/1000,fltnd); axis xy
if ~isempty(onoffset)
    nonsylzone = [[0 onoffsetm(1,1)]; [onoffsetm(1:end-1,2) onoffsetm(2:end,1)]; [onoffsetm(end,2) max(tvec)]];
    temp = [[onoffset(:,1)-margin onoffset(:,1)]; [onoffset(:,2) onoffset(:,2)+margin]];
    [~,sid] = sort(temp(:,1));
    margzone = temp(sid,:);
    idx = find((margzone(2:end,1)-margzone(1:end-1,2))>0);
    ons = [margzone(1,1);margzone(idx+1)];
    offs = [margzone(idx,2);margzone(end,2)];
    margzone = [ons offs];
else
    nonsylzone = [0 tvec(end)];
    margzone = [0 0];
end
set(0,'CurrentFigure',gv.gui.hfig);
set(gv.gui.hfig,'CurrentAxes',gv.gui.haxes(3));
hp1 = patch(margzone(:,[1 2 2 1])',repmat([0 0 fs fs]/2000,size(margzone,1),1)',1);
hp2 = patch(nonsylzone(:,[1 2 2 1])',repmat([0 0 fs fs]/2000,size(nonsylzone,1),1)',1);
set(hp1,'linestyle','none','facecolor','k','facealpha',0.1);
set(hp2,'linestyle','none','facecolor','k','facealpha',0.3);
hlf = line(tvec,freqtrace(:,1:3)/1000);
set(hlf,'linestyle','none','marker','.','color','b');
hlm = line(tvec(maxampidx),maxfreq/1000);
set(hlm,'linestyle','none','marker','+','color','r');
colormap(gca,flipud(gray));
caxis(fltndimrng);
set(gca,'tickdir','out'); box off;
ylabel('Frequency (kHz)');
xlabel('Time (s)');
set(gv.gui.haxes,'xlim',gv.dat.xl,'ylim',gv.dat.yl);    
drawnow;
% ////////////////////////////////////////////////////////////////////////
function makegui
global gv
% figure
hfig = figure;
set(hfig,'MenuBar','none','NumberTitle','off','ToolBar','none','Name',gv.usevsegver);
set(hfig,'KeyPressFcn','usvseg09r2_plus(get(gcf,''CurrentCharacter''));');
set(hfig,'DeleteFcn',@finishfunc);
% axes
haxes(1) = axes;
haxes(2) = axes;
haxes(3) = axes;
set(haxes,'tickdir','out')
% uicontrol
hpush_open = uicontrol(hfig,'style','pushbutton');
hpush_flat = uicontrol(hfig,'style','pushbutton');
hpush_thrs = uicontrol(hfig,'style','pushbutton');
hpush_dtct = uicontrol(hfig,'style','pushbutton');
hpush_save = uicontrol(hfig,'style','pushbutton');
hpush_sgmt = uicontrol(hfig,'style','pushbutton');
hpush_play = uicontrol(hfig,'style','pushbutton');
hpush_swav = uicontrol(hfig,'style','pushbutton');
hpush_long = uicontrol(hfig,'style','pushbutton');
hpush_fldr = uicontrol(hfig,'style','pushbutton');
hpush_fldr_usvcam = uicontrol(hfig,'style','pushbutton');
hedit_step = uicontrol(hfig,'style','edit');
hedit_fmin = uicontrol(hfig,'style','edit');
hedit_fmax = uicontrol(hfig,'style','edit');
hedit_thre = uicontrol(hfig,'style','edit');
hedit_dmin = uicontrol(hfig,'style','edit');
hedit_dmax = uicontrol(hfig,'style','edit');
hedit_gmin = uicontrol(hfig,'style','edit');
hedit_marg = uicontrol(hfig,'style','edit');
htggl_wavo = uicontrol(hfig,'style','toggle');
htggl_imgo = uicontrol(hfig,'style','toggle');
htggl_imgt = uicontrol(hfig,'style','toggle');
htggl_trac = uicontrol(hfig,'style','toggle');
hchk_overwriteresult = uicontrol(hfig,'style','checkbox');
hedit_read = uicontrol(hfig,'style','edit');
hedit_mapL = uicontrol(hfig,'style','edit');
hedit_mapH = uicontrol(hfig,'style','edit');
htext_fnam = uicontrol(hfig,'style','text');
htext_fslb = uicontrol(hfig,'style','text');
htext_fsvl = uicontrol(hfig,'style','text');
htext_step = uicontrol(hfig,'style','text');
htext_fmin = uicontrol(hfig,'style','text');
htext_fmax = uicontrol(hfig,'style','text');
htext_dmin = uicontrol(hfig,'style','text');
htext_dmax = uicontrol(hfig,'style','text');
htext_gmin = uicontrol(hfig,'style','text');
htext_thre = uicontrol(hfig,'style','text');
htext_marg = uicontrol(hfig,'style','text');
htext_wavo = uicontrol(hfig,'style','text');
htext_imgo = uicontrol(hfig,'style','text');
htext_imgt = uicontrol(hfig,'style','text');
htext_trac = uicontrol(hfig,'style','text');
htext_maps = uicontrol(hfig,'style','text');
htext_read = uicontrol(hfig,'style','text');
% busy mark
hpanl_busy = uipanel(hfig);
set(hpanl_busy,'backgroundcolor',[0.6 0.6 0.6]);
% strings
set(hpush_open,'string','open');
set(hpush_flat,'string','flatten');
set(hpush_thrs,'string','threshold');
set(hpush_dtct,'string','detect');
set(hpush_save,'string','save csv');
set(hpush_sgmt,'string','segment');
set(hpush_play,'string','play');
set(hpush_swav,'string','save');
set(hpush_long,'string','long file');
set(hpush_fldr,'string','multiple files');
set(hpush_fldr_usvcam,'string','usvcam files');
set(htext_fnam,'string','');
set(hedit_step,'string',num2str(gv.prm.timestep*1000));
set(hedit_fmin,'string',num2str(gv.prm.freqmin/1000));
set(hedit_fmax,'string',num2str(gv.prm.freqmax/1000));
set(hedit_thre,'string',num2str(gv.prm.threshval));
set(hedit_dmin,'string',num2str(gv.prm.durmin*1000));
set(hedit_dmax,'string',num2str(gv.prm.durmax*1000));
set(hedit_gmin,'string',num2str(gv.prm.gapmin*1000));
set(hedit_marg,'string',num2str(gv.prm.margin*1000));
set(hedit_read,'string',num2str(gv.prm.readsize));
set(hedit_mapL,'string',num2str(gv.prm.mapL/1000));
set(hedit_mapH,'string',num2str(gv.prm.mapH/1000));
set(htext_fslb,'string',{'sampling', '(kHz)'});
set(htext_fsvl,'string','-');
set(htext_step,'string',{'time step','(ms)'});
set(htext_fmin,'string',{'freq min','(kHz)'});
set(htext_fmax,'string',{'freq max','(kHz)'});
set(htext_thre,'string',{'threshold','(SD)'});
set(htext_dmin,'string',{'dur min','(ms)'});
set(htext_dmax,'string',{'dur max','(ms)'});
set(htext_gmin,'string',{'gap min','(ms)'});
set(htext_marg,'string',{'margin','(ms)'});
set(htext_wavo,'string',{'wavfile','output'});
set(htext_imgo,'string',{'image','output'});
set(htext_imgt,'string',{'image','type'});
set(htext_trac,'string',{'trace','output'});
set(htext_maps,'string',{'map','(kHz)'});
set(htext_read,'string',{'read size','(s)'});
set(htext_fnam,'string','File: filename');
% toggle
str1 = {'off','on'};
str2 = {'orig','flat'};
str3 = {'skip','redo'};
set(htggl_wavo,'value',gv.prm.wavfileoutput,'string',str1{gv.prm.wavfileoutput+1});
set(htggl_imgo,'value',gv.prm.imageoutput,'string',str1{gv.prm.imageoutput+1});
set(htggl_imgt,'value',gv.prm.imagetype,'string',str2{gv.prm.imagetype+1});
set(htggl_trac,'value',gv.prm.traceoutput,'string',str1{gv.prm.traceoutput+1});
set(hchk_overwriteresult,'value',gv.prm.overwriteresult,'string','overwrite');
set(htggl_wavo,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_imgo,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_imgt,'callback','str={''orig'',''flat''}; set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_trac,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
% callbacks
set(hpush_open,'callback','usvseg09r2_plus(''open'');');
set(hpush_flat,'callback','usvseg09r2_plus(''flat'');');
set(hpush_thrs,'callback','usvseg09r2_plus(''thrs'');');
set(hpush_dtct,'callback','usvseg09r2_plus(''dtct'');');
set(hpush_save,'callback','usvseg09r2_plus(''save'');');
set(hpush_sgmt,'callback','usvseg09r2_plus(''sgmt'');');
set(hpush_play,'callback','usvseg09r2_plus(''play'');');
set(hpush_swav,'callback','usvseg09r2_plus(''swav'');');
set(hpush_long,'callback','usvseg09r2_plus(''long'');');
set(hpush_fldr,'callback','usvseg09r2_plus(''fldr'');');
set(hpush_fldr_usvcam,'callback','usvseg09r2_plus(''usvcam'');');
% text set
set(htext_fnam,'horizontalalign','left');
% set position
set(hfig,'position',[10 80 1000 750]);
hts = [htext_fslb,htext_step,htext_fmin,htext_fmax,htext_thre,htext_dmin,htext_dmax,htext_gmin,htext_marg,htext_wavo,htext_imgo,htext_imgt,htext_trac,htext_read];
hed = [htext_fsvl,hedit_step,hedit_fmin,hedit_fmax,hedit_thre,hedit_dmin,hedit_dmax,hedit_gmin,hedit_marg,htggl_wavo,htggl_imgo,htggl_imgt,htggl_trac,hedit_read];
for n=1:length(hts)
    set(hts(n),'units','normalized','position',[0.02 0.94-(n-1)*0.035 0.06 0.035]);
    set(hed(n),'units','normalized','position',[0.08 0.94-(n-1)*0.035 0.04 0.035]);
end
set(hpush_open,'units','normalized','position',[0.02 0.42 0.10 0.04]);
set(hpush_flat,'units','normalized','position',[0.02 0.38 0.10 0.04]);
set(hpush_thrs,'units','normalized','position',[0.02 0.34 0.10 0.04]);
set(hpush_dtct,'units','normalized','position',[0.02 0.30 0.10 0.04]);
set(hpush_save,'units','normalized','position',[0.02 0.26 0.10 0.04]);
set(hpush_sgmt,'units','normalized','position',[0.02 0.22 0.10 0.04]);
set(htext_maps,'units','normalized','position',[0.02 0.16 0.04 0.04]);
set(hedit_mapL,'units','normalized','position',[0.06 0.16 0.03 0.04]);
set(hedit_mapH,'units','normalized','position',[0.09 0.16 0.03 0.04]);
set(hpush_play,'units','normalized','position',[0.02 0.12 0.05 0.04]);
set(hpush_swav,'units','normalized','position',[0.07 0.12 0.05 0.04]);
set(hpush_long,'units','normalized','position',[0.02 0.06 0.10 0.04]);
set(hpush_fldr,'units','normalized','position',[0.02 0.02 0.10 0.04]);
set(hpush_fldr_usvcam,'units','normalized','position',[0.02 0.02 0.10 0.04]);
set(hchk_overwriteresult,'units','normalized','position',[0.13 0.01 0.10 0.035]);
% filename and axis
set(hpanl_busy,'units','normalized','position',[0.175 0.955 0.02 0.025]);
set(htext_fnam,'units','normalized','position',[0.20 0.96 0.77 0.02]);
set(haxes(1), 'units','normalized','position',[0.20 0.68 0.77 0.25]);
set(haxes(2), 'units','normalized','position',[0.20 0.38 0.77 0.25]);
set(haxes(3), 'units','normalized','position',[0.20 0.08 0.77 0.25]);
% disable
set(hpush_open,'enable','on');
set(hpush_flat,'enable','off');
set(hpush_thrs,'enable','off');
set(hpush_dtct,'enable','off');
set(hpush_save,'enable','off');
set(hpush_sgmt,'enable','off');
set(hpush_play,'enable','off');
set(hpush_swav,'enable','off');
set(hpush_long,'enable','on');
set(hpush_fldr,'enable','off'); set(hpush_fldr,'visible','off');
set(hpush_fldr_usvcam,'enable','on');
% save handles
gv.gui.hfig = hfig;
gv.gui.haxes = haxes;
gv.gui.hpush_open = hpush_open;
gv.gui.hpush_flat = hpush_flat;
gv.gui.hpush_thrs = hpush_thrs;
gv.gui.hpush_dtct = hpush_dtct;
gv.gui.hpush_save = hpush_save;
gv.gui.hpush_sgmt = hpush_sgmt;
gv.gui.hpush_play = hpush_play;
gv.gui.hpush_swav = hpush_swav;
gv.gui.hpush_long = hpush_long;
gv.gui.hpush_fldr = hpush_fldr;
gv.gui.hpush_fldr_usvcam = hpush_fldr_usvcam;
gv.gui.htext_fnam = htext_fnam;
gv.gui.htext_fsvl = htext_fsvl;
gv.gui.hedit_step = hedit_step;
gv.gui.hedit_fmin = hedit_fmin;
gv.gui.hedit_fmax = hedit_fmax;
gv.gui.hedit_thre = hedit_thre;
gv.gui.hedit_dmin = hedit_dmin;
gv.gui.hedit_dmax = hedit_dmax;
gv.gui.hedit_gmin = hedit_gmin;
gv.gui.hedit_marg = hedit_marg;
gv.gui.hedit_read = hedit_read;
gv.gui.hedit_mapL = hedit_mapL;
gv.gui.hedit_mapH = hedit_mapH;
gv.gui.htggl_wavo = htggl_wavo;
gv.gui.htggl_imgo = htggl_imgo;
gv.gui.htggl_imgt = htggl_imgt;
gv.gui.htggl_trac = htggl_trac;
gv.gui.hchk_overwriteresult = hchk_overwriteresult;
gv.gui.hpanl_busy = hpanl_busy;
% delete function
function finishfunc(~,~)
global gv
prmname = [fileparts(mfilename('fullpath')) filesep 'usvseg_prm.mat'];
fetchparams;
prm = gv.prm;
save(prmname,'prm');
clearvars -global gv;
disp('bye!')
% ////////////////////////////////////////////////////////////////////////
function busytoggle(num)
global gv
if num==1
    set(gv.gui.hpanl_busy,'backgroundcolor',[1 0 0]);
    drawnow;
else
    set(gv.gui.hpanl_busy,'backgroundcolor',[0 1 0]);
    drawnow;
end 

% ////////////////////////////////////////////////////////////////////////
function fetchparams
global gv
prm.fftsize = gv.prm.fftsize;
prm.timestep = str2num(get(gv.gui.hedit_step,'string'))/1000;
prm.freqmin = str2num(get(gv.gui.hedit_fmin,'string'))*1000;
prm.freqmax = str2num(get(gv.gui.hedit_fmax,'string'))*1000;
prm.threshval = str2num(get(gv.gui.hedit_thre,'string'));
prm.durmin = str2num(get(gv.gui.hedit_dmin,'string'))/1000;
prm.durmax = str2num(get(gv.gui.hedit_dmax,'string'))/1000;
prm.gapmin = str2num(get(gv.gui.hedit_gmin,'string'))/1000;
prm.margin = str2num(get(gv.gui.hedit_marg,'string'))/1000;
prm.wavfileoutput= get(gv.gui.htggl_wavo,'value');
prm.imageoutput = get(gv.gui.htggl_imgo,'value');
prm.imagetype = get(gv.gui.htggl_imgt,'value');
prm.traceoutput = get(gv.gui.htggl_trac,'value');
prm.overwriteresult = get(gv.gui.hchk_overwriteresult,'value');
prm.readsize = str2num(get(gv.gui.hedit_read,'string'));
prm.mapL = str2num(get(gv.gui.hedit_mapL,'string'))*1000;
prm.mapH = str2num(get(gv.gui.hedit_mapH,'string'))*1000;
gv.prm = prm;

% ////////////////////////////////////////////////////////////////////////
function [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,params,med,thresh)
timestep = params.timestep;
gapmin = params.gapmin;
durmin = params.durmin;
durmax = params.durmax;
margin = params.margin;
freqmin = params.freqmin;
freqmax = params.freqmax;
% multitaper spec
mtsp = multitaperspec(wav,fs,fftsize,timestep,1);
% flattening
[fltnd,med] = flattening(mtsp,med);
% threshold calculation with n*sigma (SD) of background noise 
if isempty(thresh)
    thresh = estimatethresh(fltnd,fs,freqmin,freqmax,params.threshval);
end
% thresholding
thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh);
% onset/offset detection
[onoffset,onoffsetm,~,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
% peak tracking
[freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin);

% ////////////////////////////////////////////////modified by JM, 2021/4/12
function segfun(startid,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,tvec,amptrace,freqtrace,wavflg,imgflg,trcflg)
fftsize = (size(inputimg,1)-1)*2;
step = round(timestep*fs);
onset = onoffset(:,1);
offset = onoffset(:,2);
im = flipud(uint8((inputimg-imrng(1))/(imrng(2)-imrng(1))*64));
cm = flipud(gray(64));    

for n=1:length(onset)
    rng = [max(round((onset(n)-margin)*fs),1) min(round((offset(n)+margin)*fs),length(wav))];
    rngs = round((rng-fftsize/2)/step);
    rng2 = [max(rngs(1),1) min(rngs(2),size(im,2))];
    % wave write
    if wavflg==1
        fname = sprintf('%s_%04d.wav',prefix,n+startid-1);
        audiowrite([outp filesep fname],wav(rng(1):rng(2)),fs);
    end
    % jpg write
    if imgflg==1
        imseg = im(:,rng2(1):rng2(2));
        fname = sprintf('%s_%04d.jpg',prefix,n+startid-1);
        imwrite(imseg,cm,[outp filesep fname]);
    end
    % trace write
    if trcflg==1
        af = [amptrace(:,1) freqtrace(:,1) amptrace(:,2) freqtrace(:,2) amptrace(:,3) freqtrace(:,3)];
        dat = [tvec(rng2(1):rng2(2))' af(rng2(1):rng2(2),:)];
        % CSV file
        fname = sprintf('%s_%04d.ori.csv',prefix,n+startid-1);
        dlmwrite([outp filesep fname],dat,'precision',10);
    end
    
    % get & output subsegmens  %%%% added by JM, 2021/4/12
    if trcflg==1
        ft = freqtrace(rng2(1):rng2(2),:)*fftsize/fs + 1;
        at = amptrace(rng2(1):rng2(2),:);
        t = tvec(rng2(1):rng2(2));

        imseg = im(:,rng2(1):rng2(2));
        %outSubSegmentData(imseg, ft, at, [outp filesep sprintf('%s_%04d.ss.mat',prefix,n+startid-1)]);
        ss = getSubSegment(imseg, ft, at);
        ss(:,1) = t(ss(:,1));
        ss(:,2) = (ss(:,2)-1)*fs/fftsize;
        
        % mark subsegment totally within the margins 
        for i = 1:max(ss(:,4))
            I = ss(:,4) == i;
            if sum(I)==0
                continue;
            end
            if (max(ss(I,1)) < min(t)+margin) || (min(ss(I,1)) > max(t)-margin)
                ss(I,4) = -2;
            end
        end
        I = (ss(:,4) == -1) & ((ss(:,1) < min(t)+margin) | (ss(:,1) > max(t)-margin));
        ss(I,4) = -2;
        
        % exclude very weak amplitude segments
        for i = 1:max(ss(:,4))
            I = ss(:,4) == i;
            if sum(I)==0
                continue;
            end
            if median(ss(I,3)) < 3  % 3dB from baseline
                ss(I,4) = -3;
            end
        end
        
        fname = sprintf('%s_%04d.ss.csv',prefix,n+startid-1);

        dlmwrite([outp filesep fname],ss,'precision',10);
    end
    
end

% ////////////////////////////////////////modified by J Matsumoto 2021/4/12
function outSubSegmentData(fltnd, peakfreqsg, peakampsg, outfn)

img = zeros(size(fltnd));
for i=1:4
    f = peakfreqsg(:,i);
    t = [1:size(f,1)]';
    I = ~isnan(f);
    t = floor(t(I));
    f = floor(f(I));
    I = sub2ind(size(fltnd),f,t);
    img(I) = 1;
end
save(outfn, 'img', 'fltnd');


function ss = getSubSegment(fltnd, peakfreqsg, peakampsg)

img = zeros(size(fltnd));
for i=1:4
    f = peakfreqsg(:,i);
    t = [1:size(f,1)]';
    I = ~isnan(f);
    t = floor(t(I));
    f = floor(f(I));
    I = sub2ind(size(fltnd),f,t);
    img(I) = 1;
end

se = ones(3,3);
img = imdilate(img,se);

img = imdilate(img,se);
img = imerode(img,se);

C = corner(img, 'SensitivityFactor', 0.2, 'FilterCoefficients', fspecial('gaussian',[51 1],2.5) );
for i=1:size(C,1)
    img(max(1,C(i,2)-7):min(size(img,1),C(i,2)+7), max(1, C(i,1)-1):min(size(img,2), C(i,1)+1)) = 0;
end

L = watershed(-img);
B =L==0;

P = [];
for i=1:4
    f = peakfreqsg(:,i);
    a = peakampsg(:,i);
    t = [1:size(f,1)]';
    I = ~isnan(f);
    P = [P; t(I), f(I), a(I)];
end

A = zeros(size(P,1),1);
for i=1:size(P,1)
    A(i) = L(floor(P(i,2)), floor(P(i,1)));
end

P = [P A];

I = P(:,4)==0;
P(I,:) = [];

for i=1:max(P(:,4))
   I = P(:,4) == i;
   if sum(I)<7
       P(I,4) = -1;
   end
end
P = sortrows(P, 4);

ss = P;


% ////////////////////////////////////////modified by J Matsumoto 2021/4/12
function [freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin)
if isempty(onoffset)
    freqtrace = nan(size(mtsp,2),4);
    amptrace =  nan(size(mtsp,2),4);
    maxampval = [];
    maxampidx = [];
    maxfreq = [];
    meanfreq = [];
    cvfreq = [];
    return;
end
fftsize = (size(fltnd,1)-1)*2;
step = round(fs*timestep);
% spectral peak saliency
spcsal = spectralsaliency(fltnd);
% get peak and reorganize
bandwidth = 4;%9;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% modified by JM, 2021/4/12 
ncandidates = 4;
contmin = 10;
ampthr = 0;
nstep = size(fltnd,2);
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
onidx = round((onoffset(:,1)-margin)*fs/step); % add margin
offidx = round((onoffset(:,2)+margin)*fs/step); % add margin
onidx(1) = max(1,onidx(1));
offidx(end) = min(nstep,offidx(end));
freqmat = nan(nstep,ncandidates);
ampmat = nan(nstep,ncandidates);
maxampval = nan(length(onoffset(:,1)),1);
maxampidx = ones(length(onoffset(:,1)),1);
maxampfreq = nan(size(onoffset,1),1);
meanfreq = nan(size(onoffset,1),1);
cvfreq = nan(size(onoffset,1),1);

for n=1:size(onoffset,1)
    idx = onidx(n):offidx(n);
    [peakfreq,peakamp] = searchpeak(spcsal(fnmin:fnmax,idx),fltnd(fnmin:fnmax,idx),ncandidates,bandwidth);
    [peakfreqsg,peakampsg] = segregatepeak(peakfreq+fnmin-1,peakamp,contmin,ampthr);
  
    freqmat(idx,:) = peakfreqsg;
    ampmat(idx,:) = peakampsg;
    if any(~isnan(peakampsg(:)))
        [mvC,miC] = max(peakampsg,[],2);
        [~,miR] = max(mvC);
        maxampidx(n) = miR+idx(1)-1;
        maxampfreq(n) = peakfreqsg(miR,miC(miR));
        if ~isnan(maxampfreq(n))
            maxampval(n) = mtsp(round(maxampfreq(n)),maxampidx(n));
        end
    end
    meanfreq(n) = mean((freqmat(idx,1)-1)/fftsize*fs,'omitnan');
    ft = (peakfreqsg(:,1)-1)/fftsize*fs;
    cvfreq(n) = std(ft,1,'omitnan')/meanfreq(n);
end
freqtrace = (freqmat-1)/fftsize*fs;
amptrace = ampmat;
maxfreq = (maxampfreq-1)/fftsize*fs;

% ////////////////////////////////////////////////////////////////////////
function snd = soundsynthesis(freq,amp,tvec,fs,rfs,freqmap)
% time vector
rt = (1:round(rfs*max(tvec)))'/rfs;
% process frequency
nid = ~isnan(freq);
npf = freq(nid);
npf = [mean(npf); npf; mean(npf)];
nT = [1/fs; tvec(nid)'; max(tvec)];
p = interp1(nT,npf,rt);
pm = p/(fs/2)*(freqmap(2)-freqmap(1))+freqmap(1);
pm(pm<100) = 100; % lower limit
% process amplitude
a2 = amp;
a2(isnan(a2)) = -120;
a3 = interp1(tvec,a2,rt);
a4 = 10.^(a3/20);
afil = filter(ones(128,1)/128,1,[a4;zeros(64,1)]);
afil = afil(65:end);
ampli = 0.2 * afil/max(afil);
ampli(20*log10(afil)<-120) = 0;
% synthesis
omega = 2*pi*pm;
ph = cumsum(omega/rfs);
sig = sin(ph);
snd = sig.*ampli + 2^-16*randn(size(sig));

% ////////////////////////////////////////////////////////////////////////
function [peakfreq,peakamp] = searchpeak(specsaliency,specamp,ncandidates,bandwidth)
num_steps = size(specsaliency,2);
search_range = bandwidth-1;
remove_range = bandwidth*2-1;
peakfreq = nan(num_steps,ncandidates);
peakamp = nan(num_steps,ncandidates);
specsaliency(specsaliency<0) = 0;
for n=1:num_steps
    slice = specsaliency(:,n);
    for c=1:ncandidates
        [~,mi] = max(slice);
        % center of gravity
        rng = max(mi-search_range,1):min(mi+search_range,size(slice,1));
        temp = specsaliency(rng,n);
        peak = sum(temp.*rng')/sum(temp); 
        % store
        peakfreq(n,c) = peak;
        peakamp(n,c) = specamp(mi,n);
        % remove
        idx = max(round(peak)-remove_range,1):min(round(peak)+remove_range,size(slice,1));
        slice(idx) = -Inf;
    end
end
% ////////////////////////////////////////////////////////////////////////
function [peakfreqsg,peakampsg] = segregatepeak(peakfreq,peakamp,conthr,ampthr)
% amplitude thresholding
peakfreq(peakamp<ampthr) = nan;
peakamp(peakamp<ampthr) = nan;
%object segregatin with putting object number
%allow skipping two frames (steps)
distthr = 0.05; % 5 percent: fixed parameter
[nstep,ncand] = size(peakfreq);
objmat = reshape((1:(nstep*ncand)),ncand,nstep)';
objmat(isnan(peakfreq)) = nan;
nskip = 2; % can skip max 2 frames if intermediate framse are NaN
distmat = nan(nstep-3,ncand,nskip+1);
pathmat = nan(nstep-3,ncand,nskip+1);
for n=1:nstep-nskip-1
    for m=1:nskip+1
        temp = abs((1./peakfreq(n,:)'*peakfreq(n+m,:))-1);
        [mv,mid] = min(temp,[],2);
        distmat(n,:,m) = mv;
        pathmat(n,:,m) = mid;
    end
    pm = pathmat;
    pm(distmat>distthr) = nan;
    pm(isnan(distmat)) = nan;
    pp = pm(:,:,1);
    if any(~isnan(pm(n,:,1)))
        pp = pm(:,:,1);
        x = n+1;
    elseif any(~isnan(pm(n,:,2)))
        pp = pm(:,:,2);
        x = n+2;
    elseif any(~isnan(pm(n,:,3)))
        pp = pm(:,:,3);
        x = n+3;
    end
    for m=1:ncand
        if ~isnan(pp(n,m))
            if objmat(x,pp(n,m)) < objmat(n,m)
                val = objmat(x,pp(n,m));
                objmat(objmat==val) = objmat(n,m);
            else
                objmat(x,pp(n,m)) = objmat(n,m);
            end
        end
    end
end
% thresholding
objnum = unique(objmat(:));
objnum = objnum(~isnan(objnum));
peaks2 = peakfreq;
ampmat2 = peakamp;
objmat2 = objmat;
objlen = zeros(length(objnum),1);
objamp = zeros(length(objnum),1);
for n=1:length(objnum)
    idx = find(objmat==objnum(n));
    objlen(n) = length(idx);
    objamp(n) = mean(ampmat2(objmat==objnum(n)));
end
for n=1:length(objlen)
    if objlen(n)<conthr
        objlen(n) = nan;
        peaks2(objmat==objnum(n)) = nan;
        objmat2(objmat==objnum(n)) = nan;
        ampmat2(objmat==objnum(n)) = nan;
    end
end
objnum = objnum(~isnan(objlen));
objamp = objamp(~isnan(objlen));
objlen = objlen(~isnan(objlen));
% sorting
peakfreqsg = nan(size(peaks2));
peakampsg = nan(size(peakamp));
for n=1:nstep
    on = objmat2(n,:);
    oa = nan(length(on),1);
    for m=1:length(on)
        if ~isempty(find(objnum==on(m)))
            oa(m) = objamp(objnum==on(m));
        end
    end
    oa2 = oa;
    oa2(isnan(oa)) = -Inf;
    [~,sid] = sort(oa2,'descend');
    peakfreqsg(n,:) = peaks2(n,sid);
    peakampsg(n,:) = ampmat2(n,sid);
end

% ////////////////////////////////////////////////////////////////////////
function spcsal = spectralsaliency(fltnd)
fftsize = (size(fltnd,1)-1)*2;
tfsp = fftshift(sum(abs(fft(dpsstapers,fftsize)),2));
dtfsp = -diff(tfsp,2); % second-order differential
rng = fftsize/2+(-6:6);
rdtfsp = dtfsp(rng);
salfilt = (rdtfsp-mean(rdtfsp))/std(rdtfsp);
fil = filter(salfilt,1,[fltnd;zeros(6,size(fltnd,2))]);
spcsal = fil(7:end,:);

% ////////////////////////////////////////////////////////////////////////
function [onoffset,onoffsetm,onoffsig,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin,onoffthresh)
if nargin<8
    onoffthresh = 5; % optimized for multitaper spectrogram
end
fftsize = (size(thrshd,1)-1)*2;
step = round(timestep*fs);
% onset/offset detection
onoff = max(filter(ones(onoffthresh,1),1,thrshd))'>=onoffthresh;
% merge fragmented pieces
ndurmin = round(durmin*fs/step);
f = filter(ones(ndurmin,1)/ndurmin,1,[onoff;zeros(round(ndurmin/2),1)]);
monoff = f(round(ndurmin/2)+1:end)>0.5;
monoff(1) = 0; 
monoff(end) = onoff(end);
onidx = find(diff(monoff)>0)+1;
offidx = find(diff(monoff)<0)+1;
if isempty(onidx)||isempty(offidx)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
offidx(end) = min(offidx(end),length(monoff));
if ~isempty(onidx) && monoff(end)==1
    onidx = onidx(1:end-1);
end
% continuity flag: check if the end of read file is "ON"
contflg = monoff(end);
% gap thresholding
gap = (onidx(2:end)-offidx(1:end-1))*timestep;
gap(end+1) = 0;
gid = find(gap>=gapmin);
if ~isempty(gid)
    onidx = [onidx(1); onidx(gid+1)];
    offidx = [offidx(gid); offidx(end)];
else
    onidx = onidx(1);
    offidx = offidx(end);
end
% syllable duration threholding
dur = (offidx-onidx)/fs*step;
did = find(durmin<=dur & dur<=durmax);
onidx = onidx(did);
offidx = offidx(did);
tvec = ((0:(size(thrshd,2)-1))'*step+fftsize/2)/fs;
onset = tvec(onidx);
offset = tvec(offidx);
if isempty(onset)||isempty(offset)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
% margin addition
onsetm = onset-margin;
offsetm = offset+margin;
% syllables whose margins are overlapped are integrated in onoffsetm but not in onoffset
idx = find((onsetm(2:end)-offsetm(1:end-1))>0);
onsetI = [onset(1);onset(idx+1)];
offsetI = [offset(idx);offset(end)];
onsetm = onsetI-margin;
onsetm(1) = max(1/fs*step,onsetm(1));
offsetm = offsetI+margin;
offsetm(end) = min(max(size(thrshd,2)*step/fs),offsetm(end));    
% output 
onoffset = [onset offset];
onoffsetm = [onsetm offsetm];
% on/off signal
temp = zeros(size(onoff));
onidx2 = round((onset*fs-fftsize/2)/step+1);
offidx2 = round((offset*fs-fftsize/2)/step+1);
temp(onidx2) = 1;
temp(offidx2+1) = -1;
onoffsig = cumsum(temp);

% ////////////////////////////////////////////////////////////////////////
function thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
mask = zeros(size(fltnd));
mask(fnmin:fnmax,:) = 1;
thrshd = (fltnd>thresh).*mask;

% ////////////////////////////////////////////////////////////////////////
function thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
cut = fltnd(fnmin:fnmax,:);
bin = -0.05:0.1:10;
bc = bin(1:end-1)+diff(bin)/2;
h = histcounts(cut(:),bin);
fwhm = bc(find(h<h(1)/2,1)-1)*2;
sigma = fwhm/2.35;
thresh = sigma*threshval;

% ////////////////////////////////////////////////////////////////////////
function [fltnd,med] = flattening(mtsp,med)
% generate flattned spectrogram
if nargin<2
    med = [];
end
% liftering
liftercutoff = 6;%%%%%% changed by JM, May 7th, 2021
fftsize = (size(mtsp,1)-1)*2;
cep = fft([mtsp;flipud(mtsp(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);
% median centering
if isempty(med)
    med = median(liftered,2);
end
liftmed = liftered-repmat(med,1,size(liftered,2));
% 5-point median filter on frequency axis
if exist('movmedian')==2
    fltnd = movmedian(liftmed,5);
else % for R2015b or earlier
    pad1 = repmat(liftmed(1,:),2,1);
    pad2 = repmat(liftmed(end,:),2,1);
    padded = [pad1;liftmed;pad2];
    fltnd = zeros(size(liftmed));
    for n=1:size(liftmed,1)
        fltnd(n,:) = median(padded(n+(0:4),:));
    end
end

% ////////////////////////////////////////////////////////////////////////
function mtsp = multitaperspec(wav,fs,fftsize,timestep,parflag)
% generate multitaper spectrogram
step = round(timestep*fs);
wavlen = length(wav);
tapers = dpsstapers;
ntapers = size(tapers,2);
nsteps = floor((wavlen-fftsize+step)/step);
spgsize = fftsize/2+1;
idx = repmat((1:fftsize)',1,nsteps)+repmat((0:nsteps-1)*step,fftsize,1);
wavslice = wav(idx);
spmat = zeros(spgsize,nsteps);
if parflag == 1
    % use "parfor" if Parallel Computing Toolbox installed
    parfor n=1:ntapers 
        ft = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
        spmat = spmat + abs(ft(1:(fftsize/2+1),:));
    end
else
    for n=1:ntapers
        ft = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
        spmat = spmat + abs(ft(1:(fftsize/2+1),:));
    end
end
mtsp = 20*log10((spmat/ntapers)*sqrt(1/(2*pi*fftsize)));

% ////////////////////////////////////////////////////////////////////////
function D = dpsstapers
% DPSS windows for multitaper method
%  same as D = dpss(512,3,6)
D = [
7.5021e-05  0.00059175   0.0031111    0.012305     0.03677     0.07885
8.6739e-05  0.00066363   0.0033883    0.013036    0.037988    0.079687
9.9339e-05  0.00073958    0.003676    0.013781    0.039203    0.080476
0.00011286   0.0008197   0.0039744     0.01454    0.040415    0.081217
0.00012734  0.00090412   0.0042837    0.015313     0.04162    0.081907
0.00014281  0.00099295   0.0046039      0.0161     0.04282    0.082547
0.00015933   0.0010863    0.004935    0.016899    0.044012    0.083136
0.00017692   0.0011843   0.0052772    0.017711    0.045195    0.083673
0.00019564    0.001287   0.0056306    0.018535    0.046369    0.084157
0.00021553   0.0013946   0.0059951    0.019371    0.047532    0.084587
0.00023663   0.0015072   0.0063709    0.020219    0.048683    0.084963
0.00025898   0.0016249    0.006758    0.021077    0.049821    0.085285
0.00028264   0.0017477   0.0071564    0.021946    0.050945    0.085552
0.00030765   0.0018759   0.0075661    0.022825    0.052054    0.085763
0.00033404   0.0020095   0.0079872    0.023713    0.053147    0.085918
0.00036189   0.0021487   0.0084197    0.024611    0.054223    0.086018
0.00039122   0.0022934   0.0088636    0.025517     0.05528     0.08606
 0.0004221    0.002444   0.0093189    0.026431    0.056317    0.086046
0.00045457   0.0026004   0.0097854    0.027352    0.057334    0.085975
0.00048869   0.0027627    0.010263     0.02828     0.05833    0.085847
0.00052449   0.0029311    0.010753    0.029214    0.059303    0.085662
0.00056205   0.0031057    0.011253    0.030154    0.060252     0.08542
0.00060141   0.0032866    0.011765    0.031098    0.061177    0.085121
0.00064262   0.0034738    0.012287    0.032048    0.062076    0.084765
0.00068574   0.0036675    0.012821    0.033001    0.062949    0.084353
0.00073082   0.0038678    0.013366    0.033957    0.063794    0.083884
0.00077792   0.0040747    0.013922    0.034915    0.064611    0.083358
 0.0008271   0.0042883    0.014488    0.035875    0.065398    0.082777
0.00087841   0.0045088    0.015065    0.036837    0.066155    0.082141
0.00093191   0.0047362    0.015652    0.037799    0.066881    0.081449
0.00098766   0.0049706     0.01625     0.03876    0.067574    0.080703
 0.0010457   0.0052121    0.016858     0.03972    0.068235    0.079903
 0.0011061   0.0054608    0.017476    0.040679    0.068862     0.07905
  0.001169   0.0057167    0.018104    0.041635    0.069455    0.078144
 0.0012343   0.0059799    0.018742    0.042588    0.070013    0.077186
 0.0013022   0.0062505    0.019389    0.043538    0.070534    0.076177
 0.0013726   0.0065285    0.020045    0.044482    0.071019    0.075118
 0.0014458    0.006814    0.020711    0.045421    0.071466     0.07401
 0.0015216   0.0071071    0.021385    0.046354    0.071876    0.072853
 0.0016003   0.0074078    0.022068     0.04728    0.072247    0.071648
 0.0016818   0.0077162     0.02276    0.048198    0.072579    0.070397
 0.0017662   0.0080323    0.023459    0.049108    0.072871    0.069101
 0.0018535   0.0083561    0.024167    0.050008    0.073123     0.06776
 0.0019439   0.0086878    0.024882    0.050899    0.073334    0.066376
 0.0020374   0.0090273    0.025604    0.051779    0.073504    0.064951
 0.0021341   0.0093746    0.026333    0.052647    0.073632    0.063485
  0.002234   0.0097299    0.027069    0.053503    0.073718    0.061979
 0.0023371    0.010093    0.027811    0.054346    0.073762    0.060436
 0.0024436    0.010464     0.02856    0.055175    0.073764    0.058856
 0.0025535    0.010843    0.029314    0.055989    0.073722     0.05724
 0.0026669     0.01123    0.030073    0.056788    0.073638    0.055592
 0.0027838    0.011625    0.030837    0.057571     0.07351     0.05391
 0.0029043    0.012028    0.031606    0.058337    0.073338    0.052199
 0.0030284    0.012439    0.032379    0.059086    0.073123    0.050458
 0.0031563    0.012858    0.033156    0.059816    0.072865    0.048689
 0.0032879    0.013285    0.033937    0.060527    0.072562    0.046895
 0.0034233    0.013719    0.034721    0.061219    0.072216    0.045076
 0.0035627    0.014162    0.035507    0.061889    0.071826    0.043235
  0.003706    0.014612    0.036295    0.062539    0.071393    0.041373
 0.0038533    0.015071    0.037085    0.063167    0.070916    0.039491
 0.0040047    0.015537    0.037877    0.063772    0.070396    0.037593
 0.0041602     0.01601     0.03867    0.064354    0.069832    0.035678
 0.0043199    0.016492    0.039463    0.064912    0.069226    0.033749
 0.0044839    0.016981    0.040256    0.065445    0.068577    0.031809
 0.0046521    0.017477    0.041048    0.065954    0.067886    0.029858
 0.0048247    0.017981     0.04184    0.066436    0.067153    0.027898
 0.0050018    0.018493     0.04263    0.066892    0.066379    0.025932
 0.0051832    0.019012    0.043419    0.067321    0.065563    0.023961
 0.0053692    0.019537    0.044205    0.067723    0.064707    0.021986
 0.0055598    0.020071    0.044988    0.068096    0.063811    0.020011
  0.005755    0.020611    0.045768    0.068441    0.062875    0.018036
 0.0059548    0.021158    0.046544    0.068756      0.0619    0.016064
 0.0061594    0.021712    0.047315    0.069042    0.060887    0.014095
 0.0063687    0.022272    0.048082    0.069298    0.059837    0.012133
 0.0065828    0.022839    0.048843    0.069523    0.058749    0.010179
 0.0068018    0.023413    0.049599    0.069717    0.057626   0.0082344
 0.0070257    0.023992    0.050348    0.069879    0.056467   0.0063014
 0.0072546    0.024578     0.05109     0.07001    0.055273   0.0043817
 0.0074884     0.02517    0.051824    0.070108    0.054046   0.0024771
 0.0077273    0.025768    0.052551    0.070174    0.052786   0.0005893
 0.0079712    0.026372    0.053269    0.070206    0.051494  -0.0012799
 0.0082202    0.026981    0.053978    0.070206    0.050171  -0.0031288
 0.0084744    0.027595    0.054677    0.070172    0.048818  -0.0049556
 0.0087337    0.028215    0.055366    0.070104    0.047436  -0.0067588
 0.0089983    0.028839    0.056045    0.070002    0.046026  -0.0085367
 0.0092681    0.029469    0.056712    0.069866    0.044589   -0.010288
 0.0095431    0.030103    0.057368    0.069695    0.043127    -0.01201
 0.0098235    0.030741    0.058011     0.06949     0.04164   -0.013702
  0.010109    0.031384    0.058642     0.06925    0.040129   -0.015363
    0.0104     0.03203    0.059259    0.068976    0.038596    -0.01699
  0.010697     0.03268    0.059862    0.068666    0.037043   -0.018583
  0.010998    0.033334    0.060451    0.068321    0.035469    -0.02014
  0.011305    0.033991    0.061026    0.067942    0.033877   -0.021659
  0.011618    0.034652    0.061585    0.067527    0.032267    -0.02314
  0.011936    0.035314    0.062128    0.067078    0.030641    -0.02458
   0.01226     0.03598    0.062655    0.066593    0.029001   -0.025979
  0.012589    0.036648    0.063165    0.066074    0.027347   -0.027336
  0.012923    0.037318    0.063657     0.06552    0.025681   -0.028649
  0.013263     0.03799    0.064132    0.064932    0.024004   -0.029917
  0.013608    0.038663    0.064588    0.064309    0.022318   -0.031139
  0.013959    0.039337    0.065026    0.063652    0.020624   -0.032314
  0.014315    0.040013    0.065444    0.062961    0.018924   -0.033441
  0.014677    0.040689    0.065843    0.062237    0.017218   -0.034519
  0.015044    0.041366    0.066222    0.061479    0.015509   -0.035548
  0.015417    0.042043     0.06658    0.060688    0.013797   -0.036526
  0.015795    0.042719    0.066917    0.059864    0.012084   -0.037453
  0.016179    0.043395    0.067232    0.059008    0.010372   -0.038327
  0.016568    0.044071    0.067526    0.058119   0.0086619   -0.039149
  0.016963    0.044745    0.067798      0.0572    0.006955   -0.039918
  0.017363    0.045418    0.068047    0.056249    0.005253   -0.040633
  0.017768    0.046089    0.068273    0.055267   0.0035571   -0.041293
  0.018179    0.046758    0.068475    0.054256    0.001869   -0.041899
  0.018595    0.047425    0.068654    0.053215   0.0001899    -0.04245
  0.019016     0.04809    0.068809    0.052144  -0.0014787   -0.042945
  0.019442    0.048751     0.06894    0.051046  -0.0031353   -0.043384
  0.019874    0.049409    0.069045    0.049919  -0.0047787   -0.043768
  0.020311    0.050063    0.069126    0.048765  -0.0064073   -0.044096
  0.020753    0.050713    0.069182    0.047585  -0.0080198   -0.044368
    0.0212    0.051359    0.069211    0.046379  -0.0096149   -0.044583
  0.021652    0.052001    0.069216    0.045148   -0.011191   -0.044743
   0.02211    0.052637    0.069194    0.043892   -0.012747   -0.044847
  0.022572    0.053268    0.069145    0.042613   -0.014282   -0.044895
  0.023039    0.053894    0.069071    0.041311   -0.015794   -0.044887
  0.023511    0.054513    0.068969    0.039987   -0.017281   -0.044824
  0.023988    0.055126    0.068841    0.038642   -0.018744   -0.044707
  0.024469    0.055732    0.068685    0.037276    -0.02018   -0.044535
  0.024956    0.056331    0.068502    0.035891   -0.021588   -0.044309
  0.025446    0.056923    0.068292    0.034487   -0.022967   -0.044029
  0.025942    0.057507    0.068055    0.033065   -0.024315   -0.043697
  0.026442    0.058083    0.067789    0.031626   -0.025633   -0.043312
  0.026946    0.058651    0.067497    0.030172   -0.026918   -0.042876
  0.027454     0.05921    0.067176    0.028703    -0.02817   -0.042389
  0.027967    0.059759    0.066828     0.02722   -0.029386   -0.041852
  0.028484    0.060299    0.066451    0.025723   -0.030568   -0.041266
  0.029005     0.06083    0.066047    0.024215   -0.031712   -0.040632
   0.02953     0.06135    0.065615    0.022696    -0.03282    -0.03995
  0.030059    0.061859    0.065155    0.021167   -0.033888   -0.039222
  0.030592    0.062358    0.064668     0.01963   -0.034917   -0.038449
  0.031128    0.062846    0.064152    0.018084   -0.035906   -0.037631
  0.031668    0.063322    0.063609    0.016532   -0.036854   -0.036771
  0.032212    0.063786    0.063039    0.014974    -0.03776   -0.035869
  0.032759    0.064238    0.062441    0.013411   -0.038624   -0.034926
  0.033309    0.064678    0.061816    0.011845   -0.039444   -0.033944
  0.033862    0.065104    0.061163    0.010276    -0.04022   -0.032924
  0.034419    0.065518    0.060484   0.0087067   -0.040952   -0.031868
  0.034978    0.065918    0.059777   0.0071368   -0.041638   -0.030776
   0.03554    0.066304    0.059044   0.0055678   -0.042279    -0.02965
  0.036106    0.066676    0.058285   0.0040008   -0.042874   -0.028492
  0.036673    0.067034    0.057499   0.0024371   -0.043422   -0.027304
  0.037243    0.067376    0.056688  0.00087766   -0.043923   -0.026086
  0.037816    0.067704    0.055851 -0.00067636   -0.044376    -0.02484
  0.038391    0.068017    0.054989  -0.0022238   -0.044782   -0.023568
  0.038968    0.068314    0.054101  -0.0037636    -0.04514   -0.022272
  0.039547    0.068595    0.053189  -0.0052946    -0.04545   -0.020952
  0.040127     0.06886    0.052253  -0.0068157   -0.045711   -0.019612
   0.04071    0.069108    0.051293  -0.0083258   -0.045924   -0.018252
  0.041294     0.06934    0.050309  -0.0098238   -0.046088   -0.016874
  0.041879    0.069555    0.049302   -0.011308   -0.046203    -0.01548
  0.042466    0.069752    0.048272   -0.012779    -0.04627   -0.014071
  0.043054    0.069932     0.04722   -0.014234   -0.046289    -0.01265
  0.043643    0.070095    0.046146   -0.015672   -0.046258   -0.011218
  0.044233    0.070239    0.045051   -0.017093    -0.04618  -0.0097761
  0.044824    0.070365    0.043935   -0.018496   -0.046053  -0.0083273
  0.045415    0.070473    0.042799   -0.019878   -0.045879  -0.0068728
  0.046006    0.070562    0.041643    -0.02124   -0.045656  -0.0054144
  0.046598    0.070633    0.040467   -0.022581   -0.045387  -0.0039539
   0.04719    0.070684    0.039273   -0.023898    -0.04507   -0.002493
  0.047782    0.070716    0.038061   -0.025192   -0.044708  -0.0010336
  0.048374    0.070728    0.036831   -0.026462   -0.044299  0.00042269
  0.048965    0.070722    0.035584   -0.027705   -0.043844    0.001874
  0.049556    0.070695     0.03432   -0.028922   -0.043345   0.0033187
  0.050147    0.070648    0.033041   -0.030112   -0.042801    0.004755
  0.050736    0.070582    0.031747   -0.031273   -0.042214   0.0061811
  0.051324    0.070495    0.030438   -0.032405   -0.041583   0.0075955
  0.051912    0.070387    0.029116   -0.033507   -0.040911   0.0089963
  0.052498     0.07026     0.02778   -0.034578   -0.040196    0.010382
  0.053082    0.070112    0.026432   -0.035617   -0.039441    0.011751
  0.053665    0.069943    0.025073   -0.036623   -0.038647    0.013101
  0.054246    0.069753    0.023702   -0.037596   -0.037813    0.014432
  0.054825    0.069542    0.022321   -0.038536   -0.036941    0.015741
  0.055402    0.069311    0.020931    -0.03944   -0.036032    0.017026
  0.055977    0.069058    0.019531   -0.040309   -0.035087    0.018287
  0.056549    0.068784    0.018124   -0.041142   -0.034107    0.019522
  0.057118     0.06849    0.016709   -0.041938   -0.033093    0.020729
  0.057685    0.068174    0.015288   -0.042697   -0.032046    0.021908
  0.058249    0.067837    0.013861   -0.043418   -0.030967    0.023056
  0.058809    0.067478    0.012429     -0.0441   -0.029858    0.024172
  0.059367    0.067099    0.010993   -0.044743   -0.028719    0.025256
   0.05992    0.066698   0.0095538   -0.045347   -0.027553    0.026305
   0.06047    0.066276   0.0081118   -0.045911   -0.026359    0.027318
  0.061017    0.065833   0.0066681   -0.046435   -0.025139    0.028295
  0.061559    0.065369   0.0052234   -0.046918   -0.023896    0.029235
  0.062097    0.064883   0.0037785    -0.04736   -0.022629    0.030135
  0.062631    0.064377   0.0023343    -0.04776    -0.02134    0.030996
   0.06316    0.063849  0.00089167   -0.048119   -0.020032    0.031816
  0.063685    0.063301 -0.00054868   -0.048436   -0.018704    0.032594
  0.064205    0.062731  -0.0019859    -0.04871   -0.017359    0.033329
  0.064719    0.062141  -0.0034191   -0.048942   -0.015998    0.034022
  0.065229    0.061531  -0.0048475   -0.049132   -0.014622    0.034669
  0.065733      0.0609  -0.0062703   -0.049279   -0.013233    0.035272
  0.066232    0.060248  -0.0076866   -0.049383   -0.011832     0.03583
  0.066726    0.059576  -0.0090957   -0.049445   -0.010421    0.036341
  0.067213    0.058884   -0.010497   -0.049464  -0.0090009    0.036805
  0.067695    0.058173   -0.011889   -0.049439  -0.0075737    0.037223
   0.06817    0.057441   -0.013271   -0.049372  -0.0061407    0.037592
  0.068639     0.05669   -0.014643   -0.049263  -0.0047035    0.037914
  0.069102    0.055919   -0.016003   -0.049111  -0.0032635    0.038187
  0.069558    0.055129   -0.017351   -0.048916  -0.0018222    0.038411
  0.070008     0.05432   -0.018686    -0.04868 -0.00038114    0.038587
  0.070451    0.053492   -0.020008   -0.048401   0.0010582    0.038713
  0.070886    0.052646   -0.021314   -0.048081   0.0024942    0.038791
  0.071315    0.051781   -0.022606   -0.047719   0.0039255    0.038819
  0.071736    0.050898   -0.023881   -0.047316   0.0053505    0.038798
   0.07215    0.049998    -0.02514   -0.046872   0.0067678    0.038729
  0.072556     0.04908   -0.026381   -0.046388   0.0081759     0.03861
  0.072955    0.048144   -0.027603   -0.045865   0.0095732    0.038443
  0.073345    0.047192   -0.028806   -0.045301    0.010958    0.038227
  0.073728    0.046222    -0.02999   -0.044699     0.01233    0.037963
  0.074103    0.045237   -0.031153   -0.044059    0.013687    0.037652
  0.074469    0.044235   -0.032294    -0.04338    0.015027    0.037293
  0.074827    0.043217   -0.033414   -0.042664    0.016349    0.036888
  0.075177    0.042184   -0.034511   -0.041912    0.017652    0.036437
  0.075517    0.041136   -0.035584   -0.041124    0.018935    0.035939
   0.07585    0.040073   -0.036634     -0.0403    0.020196    0.035397
  0.076173    0.038996   -0.037659   -0.039442    0.021433    0.034811
  0.076488    0.037904   -0.038658    -0.03855    0.022646    0.034181
  0.076793    0.036799   -0.039632   -0.037625    0.023833    0.033509
  0.077089    0.035681    -0.04058   -0.036668    0.024993    0.032794
  0.077376    0.034549     -0.0415   -0.035679    0.026125    0.032039
  0.077654    0.033405   -0.042392   -0.034659    0.027228    0.031244
  0.077922    0.032249   -0.043257    -0.03361      0.0283     0.03041
  0.078181    0.031081   -0.044093   -0.032532     0.02934    0.029538
   0.07843    0.029902   -0.044899   -0.031427    0.030348     0.02863
   0.07867    0.028711   -0.045676   -0.030294    0.031322    0.027686
  0.078899     0.02751   -0.046422   -0.029135    0.032261    0.026707
  0.079119      0.0263   -0.047138   -0.027952    0.033164    0.025696
  0.079329    0.025079   -0.047823   -0.026744    0.034031    0.024652
  0.079528    0.023849   -0.048476   -0.025513     0.03486    0.023578
  0.079718    0.022611   -0.049098   -0.024261     0.03565    0.022474
  0.079898    0.021364   -0.049686   -0.022988    0.036402    0.021343
  0.080067    0.020109   -0.050243   -0.021695    0.037113    0.020185
  0.080226    0.018847   -0.050766   -0.020384    0.037783    0.019003
  0.080374    0.017577   -0.051256   -0.019055    0.038412    0.017796
  0.080513    0.016302   -0.051712    -0.01771    0.038999    0.016568
   0.08064     0.01502   -0.052134   -0.016349    0.039543    0.015319
  0.080758    0.013733   -0.052521   -0.014975    0.040044    0.014051
  0.080864     0.01244   -0.052875   -0.013588    0.040501    0.012766
  0.080961    0.011143   -0.053194   -0.012189    0.040913    0.011465
  0.081046   0.0098419   -0.053477   -0.010779    0.041281     0.01015
  0.081121   0.0085371   -0.053726  -0.0093609    0.041603    0.008822
  0.081185   0.0072291   -0.053939  -0.0079344    0.041881   0.0074833
  0.081239   0.0059184   -0.054117   -0.006501    0.042112   0.0061354
  0.081282   0.0046055    -0.05426  -0.0050621    0.042297   0.0047798
  0.081314   0.0032908   -0.054367  -0.0036189    0.042437   0.0034184
  0.081336    0.001975   -0.054438  -0.0021726     0.04253   0.0020527
  0.081346  0.00065841   -0.054474  -0.0007244    0.042576  0.00068453
  0.081346 -0.00065841   -0.054474   0.0007244    0.042576 -0.00068453
  0.081336   -0.001975   -0.054438   0.0021726     0.04253  -0.0020527
  0.081314  -0.0032908   -0.054367   0.0036189    0.042437  -0.0034184
  0.081282  -0.0046055    -0.05426   0.0050621    0.042297  -0.0047798
  0.081239  -0.0059184   -0.054117    0.006501    0.042112  -0.0061354
  0.081185  -0.0072291   -0.053939   0.0079344    0.041881  -0.0074833
  0.081121  -0.0085371   -0.053726   0.0093609    0.041603   -0.008822
  0.081046  -0.0098419   -0.053477    0.010779    0.041281    -0.01015
  0.080961   -0.011143   -0.053194    0.012189    0.040913   -0.011465
  0.080864    -0.01244   -0.052875    0.013588    0.040501   -0.012766
  0.080758   -0.013733   -0.052521    0.014975    0.040044   -0.014051
   0.08064    -0.01502   -0.052134    0.016349    0.039543   -0.015319
  0.080513   -0.016302   -0.051712     0.01771    0.038999   -0.016568
  0.080374   -0.017577   -0.051256    0.019055    0.038412   -0.017796
  0.080226   -0.018847   -0.050766    0.020384    0.037783   -0.019003
  0.080067   -0.020109   -0.050243    0.021695    0.037113   -0.020185
  0.079898   -0.021364   -0.049686    0.022988    0.036402   -0.021343
  0.079718   -0.022611   -0.049098    0.024261     0.03565   -0.022474
  0.079528   -0.023849   -0.048476    0.025513     0.03486   -0.023578
  0.079329   -0.025079   -0.047823    0.026744    0.034031   -0.024652
  0.079119     -0.0263   -0.047138    0.027952    0.033164   -0.025696
  0.078899    -0.02751   -0.046422    0.029135    0.032261   -0.026707
   0.07867   -0.028711   -0.045676    0.030294    0.031322   -0.027686
   0.07843   -0.029902   -0.044899    0.031427    0.030348    -0.02863
  0.078181   -0.031081   -0.044093    0.032532     0.02934   -0.029538
  0.077922   -0.032249   -0.043257     0.03361      0.0283    -0.03041
  0.077654   -0.033405   -0.042392    0.034659    0.027228   -0.031244
  0.077376   -0.034549     -0.0415    0.035679    0.026125   -0.032039
  0.077089   -0.035681    -0.04058    0.036668    0.024993   -0.032794
  0.076793   -0.036799   -0.039632    0.037625    0.023833   -0.033509
  0.076488   -0.037904   -0.038658     0.03855    0.022646   -0.034181
  0.076173   -0.038996   -0.037659    0.039442    0.021433   -0.034811
   0.07585   -0.040073   -0.036634      0.0403    0.020196   -0.035397
  0.075517   -0.041136   -0.035584    0.041124    0.018935   -0.035939
  0.075177   -0.042184   -0.034511    0.041912    0.017652   -0.036437
  0.074827   -0.043217   -0.033414    0.042664    0.016349   -0.036888
  0.074469   -0.044235   -0.032294     0.04338    0.015027   -0.037293
  0.074103   -0.045237   -0.031153    0.044059    0.013687   -0.037652
  0.073728   -0.046222    -0.02999    0.044699     0.01233   -0.037963
  0.073345   -0.047192   -0.028806    0.045301    0.010958   -0.038227
  0.072955   -0.048144   -0.027603    0.045865   0.0095732   -0.038443
  0.072556    -0.04908   -0.026381    0.046388   0.0081759    -0.03861
   0.07215   -0.049998    -0.02514    0.046872   0.0067678   -0.038729
  0.071736   -0.050898   -0.023881    0.047316   0.0053505   -0.038798
  0.071315   -0.051781   -0.022606    0.047719   0.0039255   -0.038819
  0.070886   -0.052646   -0.021314    0.048081   0.0024942   -0.038791
  0.070451   -0.053492   -0.020008    0.048401   0.0010582   -0.038713
  0.070008    -0.05432   -0.018686     0.04868 -0.00038114   -0.038587
  0.069558   -0.055129   -0.017351    0.048916  -0.0018222   -0.038411
  0.069102   -0.055919   -0.016003    0.049111  -0.0032635   -0.038187
  0.068639    -0.05669   -0.014643    0.049263  -0.0047035   -0.037914
   0.06817   -0.057441   -0.013271    0.049372  -0.0061407   -0.037592
  0.067695   -0.058173   -0.011889    0.049439  -0.0075737   -0.037223
  0.067213   -0.058884   -0.010497    0.049464  -0.0090009   -0.036805
  0.066726   -0.059576  -0.0090957    0.049445   -0.010421   -0.036341
  0.066232   -0.060248  -0.0076866    0.049383   -0.011832    -0.03583
  0.065733     -0.0609  -0.0062703    0.049279   -0.013233   -0.035272
  0.065229   -0.061531  -0.0048475    0.049132   -0.014622   -0.034669
  0.064719   -0.062141  -0.0034191    0.048942   -0.015998   -0.034022
  0.064205   -0.062731  -0.0019859     0.04871   -0.017359   -0.033329
  0.063685   -0.063301 -0.00054868    0.048436   -0.018704   -0.032594
   0.06316   -0.063849  0.00089167    0.048119   -0.020032   -0.031816
  0.062631   -0.064377   0.0023343     0.04776    -0.02134   -0.030996
  0.062097   -0.064883   0.0037785     0.04736   -0.022629   -0.030135
  0.061559   -0.065369   0.0052234    0.046918   -0.023896   -0.029235
  0.061017   -0.065833   0.0066681    0.046435   -0.025139   -0.028295
   0.06047   -0.066276   0.0081118    0.045911   -0.026359   -0.027318
   0.05992   -0.066698   0.0095538    0.045347   -0.027553   -0.026305
  0.059367   -0.067099    0.010993    0.044743   -0.028719   -0.025256
  0.058809   -0.067478    0.012429      0.0441   -0.029858   -0.024172
  0.058249   -0.067837    0.013861    0.043418   -0.030967   -0.023056
  0.057685   -0.068174    0.015288    0.042697   -0.032046   -0.021908
  0.057118    -0.06849    0.016709    0.041938   -0.033093   -0.020729
  0.056549   -0.068784    0.018124    0.041142   -0.034107   -0.019522
  0.055977   -0.069058    0.019531    0.040309   -0.035087   -0.018287
  0.055402   -0.069311    0.020931     0.03944   -0.036032   -0.017026
  0.054825   -0.069542    0.022321    0.038536   -0.036941   -0.015741
  0.054246   -0.069753    0.023702    0.037596   -0.037813   -0.014432
  0.053665   -0.069943    0.025073    0.036623   -0.038647   -0.013101
  0.053082   -0.070112    0.026432    0.035617   -0.039441   -0.011751
  0.052498    -0.07026     0.02778    0.034578   -0.040196   -0.010382
  0.051912   -0.070387    0.029116    0.033507   -0.040911  -0.0089963
  0.051324   -0.070495    0.030438    0.032405   -0.041583  -0.0075955
  0.050736   -0.070582    0.031747    0.031273   -0.042214  -0.0061811
  0.050147   -0.070648    0.033041    0.030112   -0.042801   -0.004755
  0.049556   -0.070695     0.03432    0.028922   -0.043345  -0.0033187
  0.048965   -0.070722    0.035584    0.027705   -0.043844   -0.001874
  0.048374   -0.070728    0.036831    0.026462   -0.044299 -0.00042269
  0.047782   -0.070716    0.038061    0.025192   -0.044708   0.0010336
   0.04719   -0.070684    0.039273    0.023898    -0.04507    0.002493
  0.046598   -0.070633    0.040467    0.022581   -0.045387   0.0039539
  0.046006   -0.070562    0.041643     0.02124   -0.045656   0.0054144
  0.045415   -0.070473    0.042799    0.019878   -0.045879   0.0068728
  0.044824   -0.070365    0.043935    0.018496   -0.046053   0.0083273
  0.044233   -0.070239    0.045051    0.017093    -0.04618   0.0097761
  0.043643   -0.070095    0.046146    0.015672   -0.046258    0.011218
  0.043054   -0.069932     0.04722    0.014234   -0.046289     0.01265
  0.042466   -0.069752    0.048272    0.012779    -0.04627    0.014071
  0.041879   -0.069555    0.049302    0.011308   -0.046203     0.01548
  0.041294    -0.06934    0.050309   0.0098238   -0.046088    0.016874
   0.04071   -0.069108    0.051293   0.0083258   -0.045924    0.018252
  0.040127    -0.06886    0.052253   0.0068157   -0.045711    0.019612
  0.039547   -0.068595    0.053189   0.0052946    -0.04545    0.020952
  0.038968   -0.068314    0.054101   0.0037636    -0.04514    0.022272
  0.038391   -0.068017    0.054989   0.0022238   -0.044782    0.023568
  0.037816   -0.067704    0.055851  0.00067636   -0.044376     0.02484
  0.037243   -0.067376    0.056688 -0.00087766   -0.043923    0.026086
  0.036673   -0.067034    0.057499  -0.0024371   -0.043422    0.027304
  0.036106   -0.066676    0.058285  -0.0040008   -0.042874    0.028492
   0.03554   -0.066304    0.059044  -0.0055678   -0.042279     0.02965
  0.034978   -0.065918    0.059777  -0.0071368   -0.041638    0.030776
  0.034419   -0.065518    0.060484  -0.0087067   -0.040952    0.031868
  0.033862   -0.065104    0.061163   -0.010276    -0.04022    0.032924
  0.033309   -0.064678    0.061816   -0.011845   -0.039444    0.033944
  0.032759   -0.064238    0.062441   -0.013411   -0.038624    0.034926
  0.032212   -0.063786    0.063039   -0.014974    -0.03776    0.035869
  0.031668   -0.063322    0.063609   -0.016532   -0.036854    0.036771
  0.031128   -0.062846    0.064152   -0.018084   -0.035906    0.037631
  0.030592   -0.062358    0.064668    -0.01963   -0.034917    0.038449
  0.030059   -0.061859    0.065155   -0.021167   -0.033888    0.039222
   0.02953    -0.06135    0.065615   -0.022696    -0.03282     0.03995
  0.029005    -0.06083    0.066047   -0.024215   -0.031712    0.040632
  0.028484   -0.060299    0.066451   -0.025723   -0.030568    0.041266
  0.027967   -0.059759    0.066828    -0.02722   -0.029386    0.041852
  0.027454    -0.05921    0.067176   -0.028703    -0.02817    0.042389
  0.026946   -0.058651    0.067497   -0.030172   -0.026918    0.042876
  0.026442   -0.058083    0.067789   -0.031626   -0.025633    0.043312
  0.025942   -0.057507    0.068055   -0.033065   -0.024315    0.043697
  0.025446   -0.056923    0.068292   -0.034487   -0.022967    0.044029
  0.024956   -0.056331    0.068502   -0.035891   -0.021588    0.044309
  0.024469   -0.055732    0.068685   -0.037276    -0.02018    0.044535
  0.023988   -0.055126    0.068841   -0.038642   -0.018744    0.044707
  0.023511   -0.054513    0.068969   -0.039987   -0.017281    0.044824
  0.023039   -0.053894    0.069071   -0.041311   -0.015794    0.044887
  0.022572   -0.053268    0.069145   -0.042613   -0.014282    0.044895
   0.02211   -0.052637    0.069194   -0.043892   -0.012747    0.044847
  0.021652   -0.052001    0.069216   -0.045148   -0.011191    0.044743
    0.0212   -0.051359    0.069211   -0.046379  -0.0096149    0.044583
  0.020753   -0.050713    0.069182   -0.047585  -0.0080198    0.044368
  0.020311   -0.050063    0.069126   -0.048765  -0.0064073    0.044096
  0.019874   -0.049409    0.069045   -0.049919  -0.0047787    0.043768
  0.019442   -0.048751     0.06894   -0.051046  -0.0031353    0.043384
  0.019016    -0.04809    0.068809   -0.052144  -0.0014787    0.042945
  0.018595   -0.047425    0.068654   -0.053215   0.0001899     0.04245
  0.018179   -0.046758    0.068475   -0.054256    0.001869    0.041899
  0.017768   -0.046089    0.068273   -0.055267   0.0035571    0.041293
  0.017363   -0.045418    0.068047   -0.056249    0.005253    0.040633
  0.016963   -0.044745    0.067798     -0.0572    0.006955    0.039918
  0.016568   -0.044071    0.067526   -0.058119   0.0086619    0.039149
  0.016179   -0.043395    0.067232   -0.059008    0.010372    0.038327
  0.015795   -0.042719    0.066917   -0.059864    0.012084    0.037453
  0.015417   -0.042043     0.06658   -0.060688    0.013797    0.036526
  0.015044   -0.041366    0.066222   -0.061479    0.015509    0.035548
  0.014677   -0.040689    0.065843   -0.062237    0.017218    0.034519
  0.014315   -0.040013    0.065444   -0.062961    0.018924    0.033441
  0.013959   -0.039337    0.065026   -0.063652    0.020624    0.032314
  0.013608   -0.038663    0.064588   -0.064309    0.022318    0.031139
  0.013263    -0.03799    0.064132   -0.064932    0.024004    0.029917
  0.012923   -0.037318    0.063657    -0.06552    0.025681    0.028649
  0.012589   -0.036648    0.063165   -0.066074    0.027347    0.027336
   0.01226    -0.03598    0.062655   -0.066593    0.029001    0.025979
  0.011936   -0.035314    0.062128   -0.067078    0.030641     0.02458
  0.011618   -0.034652    0.061585   -0.067527    0.032267     0.02314
  0.011305   -0.033991    0.061026   -0.067942    0.033877    0.021659
  0.010998   -0.033334    0.060451   -0.068321    0.035469     0.02014
  0.010697    -0.03268    0.059862   -0.068666    0.037043    0.018583
    0.0104    -0.03203    0.059259   -0.068976    0.038596     0.01699
  0.010109   -0.031384    0.058642    -0.06925    0.040129    0.015363
 0.0098235   -0.030741    0.058011    -0.06949     0.04164    0.013702
 0.0095431   -0.030103    0.057368   -0.069695    0.043127     0.01201
 0.0092681   -0.029469    0.056712   -0.069866    0.044589    0.010288
 0.0089983   -0.028839    0.056045   -0.070002    0.046026   0.0085367
 0.0087337   -0.028215    0.055366   -0.070104    0.047436   0.0067588
 0.0084744   -0.027595    0.054677   -0.070172    0.048818   0.0049556
 0.0082202   -0.026981    0.053978   -0.070206    0.050171   0.0031288
 0.0079712   -0.026372    0.053269   -0.070206    0.051494   0.0012799
 0.0077273   -0.025768    0.052551   -0.070174    0.052786  -0.0005893
 0.0074884    -0.02517    0.051824   -0.070108    0.054046  -0.0024771
 0.0072546   -0.024578     0.05109    -0.07001    0.055273  -0.0043817
 0.0070257   -0.023992    0.050348   -0.069879    0.056467  -0.0063014
 0.0068018   -0.023413    0.049599   -0.069717    0.057626  -0.0082344
 0.0065828   -0.022839    0.048843   -0.069523    0.058749   -0.010179
 0.0063687   -0.022272    0.048082   -0.069298    0.059837   -0.012133
 0.0061594   -0.021712    0.047315   -0.069042    0.060887   -0.014095
 0.0059548   -0.021158    0.046544   -0.068756      0.0619   -0.016064
  0.005755   -0.020611    0.045768   -0.068441    0.062875   -0.018036
 0.0055598   -0.020071    0.044988   -0.068096    0.063811   -0.020011
 0.0053692   -0.019537    0.044205   -0.067723    0.064707   -0.021986
 0.0051832   -0.019012    0.043419   -0.067321    0.065563   -0.023961
 0.0050018   -0.018493     0.04263   -0.066892    0.066379   -0.025932
 0.0048247   -0.017981     0.04184   -0.066436    0.067153   -0.027898
 0.0046521   -0.017477    0.041048   -0.065954    0.067886   -0.029858
 0.0044839   -0.016981    0.040256   -0.065445    0.068577   -0.031809
 0.0043199   -0.016492    0.039463   -0.064912    0.069226   -0.033749
 0.0041602    -0.01601     0.03867   -0.064354    0.069832   -0.035678
 0.0040047   -0.015537    0.037877   -0.063772    0.070396   -0.037593
 0.0038533   -0.015071    0.037085   -0.063167    0.070916   -0.039491
  0.003706   -0.014612    0.036295   -0.062539    0.071393   -0.041373
 0.0035627   -0.014162    0.035507   -0.061889    0.071826   -0.043235
 0.0034233   -0.013719    0.034721   -0.061219    0.072216   -0.045076
 0.0032879   -0.013285    0.033937   -0.060527    0.072562   -0.046895
 0.0031563   -0.012858    0.033156   -0.059816    0.072865   -0.048689
 0.0030284   -0.012439    0.032379   -0.059086    0.073123   -0.050458
 0.0029043   -0.012028    0.031606   -0.058337    0.073338   -0.052199
 0.0027838   -0.011625    0.030837   -0.057571     0.07351    -0.05391
 0.0026669    -0.01123    0.030073   -0.056788    0.073638   -0.055592
 0.0025535   -0.010843    0.029314   -0.055989    0.073722    -0.05724
 0.0024436   -0.010464     0.02856   -0.055175    0.073764   -0.058856
 0.0023371   -0.010093    0.027811   -0.054346    0.073762   -0.060436
  0.002234  -0.0097299    0.027069   -0.053503    0.073718   -0.061979
 0.0021341  -0.0093746    0.026333   -0.052647    0.073632   -0.063485
 0.0020374  -0.0090273    0.025604   -0.051779    0.073504   -0.064951
 0.0019439  -0.0086878    0.024882   -0.050899    0.073334   -0.066376
 0.0018535  -0.0083561    0.024167   -0.050008    0.073123    -0.06776
 0.0017662  -0.0080323    0.023459   -0.049108    0.072871   -0.069101
 0.0016818  -0.0077162     0.02276   -0.048198    0.072579   -0.070397
 0.0016003  -0.0074078    0.022068    -0.04728    0.072247   -0.071648
 0.0015216  -0.0071071    0.021385   -0.046354    0.071876   -0.072853
 0.0014458   -0.006814    0.020711   -0.045421    0.071466    -0.07401
 0.0013726  -0.0065285    0.020045   -0.044482    0.071019   -0.075118
 0.0013022  -0.0062505    0.019389   -0.043538    0.070534   -0.076177
 0.0012343  -0.0059799    0.018742   -0.042588    0.070013   -0.077186
  0.001169  -0.0057167    0.018104   -0.041635    0.069455   -0.078144
 0.0011061  -0.0054608    0.017476   -0.040679    0.068862    -0.07905
 0.0010457  -0.0052121    0.016858    -0.03972    0.068235   -0.079903
0.00098766  -0.0049706     0.01625    -0.03876    0.067574   -0.080703
0.00093191  -0.0047362    0.015652   -0.037799    0.066881   -0.081449
0.00087841  -0.0045088    0.015065   -0.036837    0.066155   -0.082141
 0.0008271  -0.0042883    0.014488   -0.035875    0.065398   -0.082777
0.00077792  -0.0040747    0.013922   -0.034915    0.064611   -0.083358
0.00073082  -0.0038678    0.013366   -0.033957    0.063794   -0.083884
0.00068574  -0.0036675    0.012821   -0.033001    0.062949   -0.084353
0.00064262  -0.0034738    0.012287   -0.032048    0.062076   -0.084765
0.00060141  -0.0032866    0.011765   -0.031098    0.061177   -0.085121
0.00056205  -0.0031057    0.011253   -0.030154    0.060252    -0.08542
0.00052449  -0.0029311    0.010753   -0.029214    0.059303   -0.085662
0.00048869  -0.0027627    0.010263    -0.02828     0.05833   -0.085847
0.00045457  -0.0026004   0.0097854   -0.027352    0.057334   -0.085975
 0.0004221   -0.002444   0.0093189   -0.026431    0.056317   -0.086046
0.00039122  -0.0022934   0.0088636   -0.025517     0.05528    -0.08606
0.00036189  -0.0021487   0.0084197   -0.024611    0.054223   -0.086018
0.00033404  -0.0020095   0.0079872   -0.023713    0.053147   -0.085918
0.00030765  -0.0018759   0.0075661   -0.022825    0.052054   -0.085763
0.00028264  -0.0017477   0.0071564   -0.021946    0.050945   -0.085552
0.00025898  -0.0016249    0.006758   -0.021077    0.049821   -0.085285
0.00023663  -0.0015072   0.0063709   -0.020219    0.048683   -0.084963
0.00021553  -0.0013946   0.0059951   -0.019371    0.047532   -0.084587
0.00019564   -0.001287   0.0056306   -0.018535    0.046369   -0.084157
0.00017692  -0.0011843   0.0052772   -0.017711    0.045195   -0.083673
0.00015933  -0.0010863    0.004935   -0.016899    0.044012   -0.083136
0.00014281 -0.00099295   0.0046039     -0.0161     0.04282   -0.082547
0.00012734 -0.00090412   0.0042837   -0.015313     0.04162   -0.081907
0.00011286  -0.0008197   0.0039744    -0.01454    0.040415   -0.081217
9.9339e-05 -0.00073958    0.003676   -0.013781    0.039203   -0.080476
8.6739e-05 -0.00066363   0.0033883   -0.013036    0.037988   -0.079687
7.5021e-05 -0.00059175   0.0031111   -0.012305     0.03677    -0.07885
];

