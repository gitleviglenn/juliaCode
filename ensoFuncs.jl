# storehouse for functions used in the analysis of data for my
# enso and RI work.  

function fig_plot(inpv,xx,tit)
    #println(inpv)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(1000,200),
        )
    ax = Axis(f2[1,1];
        xlabel="time (monthly)",
        ylabel="temp anom (K)",
        title=tit,
        )
    #lines!(ax, 0.:1:xx-1, inpv -> cos(inpv);
    lines!(ax, 0.:1:xx-1, inpv;
    #scatterlines!(ax, xx, inpv;
        color=:black,
        linewidth=2,
        linestyle=:solid,
        )
    return f2
    #println(inpv[1:40])
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

function smooth_ts(inpTS,len)
    # smooth the time series with a running mean
    global ts_sm = zeros(len)
    istart= 2
    jend  = len
    for i in istart:jend-1
      ts_sm[i] = (inpTS[i+1]+inpTS[i]+inpTS[i-1])/3
    end
    ts_sm[1]=ts_sm[2]
    ts_sm[jend]=ts_sm[jend-1]

    return ts_sm
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


function smooth_12_ts(inpTS,len)
    # smooth the time series with a running mean
    global ts_12_sm = zeros(len)
    istart= 6
    jend  = len
    for i in istart:jend-6
      ts_12_sm[i] = (inpTS[i+6]+inpTS[i+5]+inpTS[i+4]+inpTS[i+3]+inpTS[i+2]+inpTS[i+1]+inpTS[i]+inpTS[i-1]+inpTS[i-2]+inpTS[i-3]+inpTS[i-4]+inpTS[i-5])/12
    end
    return ts_12_sm
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

function check_thresh_high(inpFile, inpTS, thshd)
    ds = NCDataset(inpFile)
    ds.attrib
    nctime = ds["time"]
    println("~~~~~~~~~Silence~~~~~~~~~~~~~~~~~~")
    println(nctime[1:10])
    println("~~~~~~~~~Golden~~~~~~~~~~~~~~~~~~")
    # for incoming file, get the time dimension 
    # for the incoming timeseries, check at each index, if the value is greater than 
    # the treshold, if so return true, otherwise return false.   
    # or perhaps return 1 or 0.
    # for i in full lenght, is val >= thshd, grab i value
    #
    global high = findall(inpTS .> thshd)
    #
    println("~~~~~~~~~high Threshold values~~~~~~~~~~~~~~~~~~")
    println(inpTS[high])
    println("~~~~~~~~~End high Threshold Values~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~high Threshold Times     ~~~~~~~~~~~~~~~~~~")
    println(nctime[high])
    println("~~~~~~~~~End high Threshold Times     ~~~~~~~~~~~~~~~~~~")
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

function check_thresh_low(inpFile, inpTS, thshd)
    ds = NCDataset(inpFile)
    ds.attrib
    nctime = ds["time"]
    # for incoming file, get the time dimension 
    # for the incoming timeseries, check at each index, if the value is greater than 
    # the treshold, if so return true, otherwise return false.   
    # or perhaps return 1 or 0.
    # for i in full lenght, is val >= thshd, grab i value
    #
    global low = findall(inpTS .< thshd)
    #
    println("~~~~~~~~~low Threshold values~~~~~~~~~~~~~~~~~~")
    println(inpTS[low])
    println("~~~~~~~~~End low Threshold Values~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~low Threshold Times     ~~~~~~~~~~~~~~~~~~")
    println(nctime[low])
    println("~~~~~~~~~low End Threshold Times     ~~~~~~~~~~~~~~~~~~")
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


function prepare_cmip_ts(inpFile,len,troplat1,troplat2,ln1,ln2,lt1,lt2)
#function prepare_cmip_ts(inpFile,len)
    # the seasonal cycle needs to be removed
    # the tropical mean in L'Heureux are defined as +/-20, it does not 
    # appear that is being taken into account for potentially different grids.  Check!
    ds1 = NCDataset(inpFile)
    ds1.attrib
    # to print meta data for a particular variable: 
    ds1["tos"]
    sst1 = ds1["tos"] # what are the dimensions here?   Should be +/-20 degrees
    nclat = ds1["lat"]
    nclon = ds1["lon"]
    nctime = ds1["time"]
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~ nino34 long points: ~~~~~~~~~~~~~~~~")
    #println(nclon[10:61])
    println(nclon[ln1:ln2])
    println("~~~~~ nino34 lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[lt1:lt2])
    println("~~~~~tropical boundary lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[troplat1:troplat2])
    #
    nino34_full = sst1[ln1:ln2, lt1:lt2, :];
    trop_mean   = sst1[:, troplat1:troplat2, :];
    #trop_mean   = sst1[:, :, :];
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #println("~~~~type of nino34_full ~~~~~~~~~~~~~~~~~~~~~~~")
    #typeof(nino34_full)
    #size(nino34_full)
    #println(nino34_full[:,7,10])
    #println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    tlength = len
    global nino34_ts = zeros(tlength)
    tropmn_ts = zeros(tlength)

    # compute climatological average value in nino34 and tropical regions 
    nino34_ts_mn = mean(filter(!isnan, skipmissing(nino34_full)))
    tropmn_ts_mn = mean(filter(!isnan, skipmissing(trop_mean)))

    # compute the anomalies
    # nino3.4 index --> nino34_ts
    # ONI           --> 3 month runnign mean of nino34_ts
    for i = 1:tlength
        #nino34_ts[i]=mean(skipmissing(nino34_full[:,:,i]))-nino34_ts_mn
        nino34_ts[i]=mean(filter(!isnan, skipmissing(nino34_full[:,:,i])))-nino34_ts_mn;
        tropmn_ts[i]=mean(filter(!isnan, skipmissing(trop_mean[:,:,i])))-tropmn_ts_mn;
    end
    # remove the seasonal cycle:
    # note: these months are not necessarily correct.  it depends on what month the ts starts on...
    ss    = zeros(12)
    sst   = zeros(12)
    jan   = [mean([nino34_ts[i] for i in 1:12:tlength])]
    janst = std(nino34_ts[i] for i in 1:12:tlength)
    feb   = [mean([nino34_ts[i] for i in 2:12:tlength])]
    febst = std(nino34_ts[i] for i in 2:12:tlength)
    mar   = [mean([nino34_ts[i] for i in 3:12:tlength])]
    marst = std(nino34_ts[i] for i in 3:12:tlength)
    apr   = [mean([nino34_ts[i] for i in 4:12:tlength])]
    aprst = std(nino34_ts[i] for i in 4:12:tlength)
    may   = [mean([nino34_ts[i] for i in 5:12:tlength])]
    mayst = std(nino34_ts[i] for i in 5:12:tlength)
    jun   = [mean([nino34_ts[i] for i in 6:12:tlength])]
    junst = std(nino34_ts[i] for i in 6:12:tlength)
    jul   = [mean([nino34_ts[i] for i in 7:12:tlength])]
    julst = std(nino34_ts[i] for i in 7:12:tlength)
    aug   = [mean([nino34_ts[i] for i in 8:12:tlength])]
    augst = std(nino34_ts[i] for i in 8:12:tlength)
    sep   = [mean([nino34_ts[i] for i in 9:12:tlength])]
    sepst = std(nino34_ts[i] for i in 9:12:tlength)
    oct   = [mean([nino34_ts[i] for i in 10:12:tlength])]
    octst = std(nino34_ts[i] for i in 10:12:tlength)
    nov   = [mean([nino34_ts[i] for i in 11:12:tlength])]
    novst = std(nino34_ts[i] for i in 11:12:tlength)
    dec   = [mean([nino34_ts[i] for i in 12:12:tlength])]
    decst = std(nino34_ts[i] for i in 12:12:tlength)

    ss    = [jan feb mar apr may jun jul aug sep oct nov dec]
    sst   = [janst febst marst aprst mayst junst julst augst sepst octst novst decst]

    # computed ts for nino 3.4 minus seasonal cycle
    global ts_rmn_nsc = zeros(tlength)
    for i in 1:12:tlength
      ts_rmn_nsc[i]    = nino34_ts[i]    - ss[1]
      ts_rmn_nsc[i+1]  = nino34_ts[i+1]  - ss[2]
      ts_rmn_nsc[i+2]  = nino34_ts[i+2]  - ss[3]
      ts_rmn_nsc[i+3]  = nino34_ts[i+3]  - ss[4]
      ts_rmn_nsc[i+4]  = nino34_ts[i+4]  - ss[5]
      ts_rmn_nsc[i+5]  = nino34_ts[i+5]  - ss[6]
      ts_rmn_nsc[i+6]  = nino34_ts[i+6]  - ss[7]
      ts_rmn_nsc[i+7]  = nino34_ts[i+7]  - ss[8]
      ts_rmn_nsc[i+8]  = nino34_ts[i+8]  - ss[9]
      ts_rmn_nsc[i+9]  = nino34_ts[i+9]  - ss[10]
      ts_rmn_nsc[i+10] = nino34_ts[i+10] - ss[11]
      ts_rmn_nsc[i+11] = nino34_ts[i+11] - ss[12]
    end
#-------------------------------------------------------------------------- 
    # calculate ONI: compute the 3mn running average
    global ts_oni_cmip        = zeros(tlength)
    istart= 2
    jend  = tlength
    for i in istart:jend-1
      ts_oni_cmip[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3;
    end
    ts_oni_cmip[1]=ts_oni_cmip[2]
    ts_oni_cmip[jend]=ts_oni_cmip[jend-1]

    # now compute the standard deviation of ONI
    st_oni    = zeros(12)
    janst_oni = std(ts_oni_cmip[i] for i in 1:12:tlength)
    febst_oni = std(ts_oni_cmip[i] for i in 2:12:tlength)
    marst_oni = std(ts_oni_cmip[i] for i in 3:12:tlength)
    aprst_oni = std(ts_oni_cmip[i] for i in 4:12:tlength)
    mayst_oni = std(ts_oni_cmip[i] for i in 5:12:tlength)
    junst_oni = std(ts_oni_cmip[i] for i in 6:12:tlength)
    julst_oni = std(ts_oni_cmip[i] for i in 7:12:tlength)
    augst_oni = std(ts_oni_cmip[i] for i in 8:12:tlength)
    sepst_oni = std(ts_oni_cmip[i] for i in 9:12:tlength)
    octst_oni = std(ts_oni_cmip[i] for i in 10:12:tlength)
    novst_oni = std(ts_oni_cmip[i] for i in 11:12:tlength)
    decst_oni = std(ts_oni_cmip[i] for i in 12:12:tlength)
    st_oni    = [janst_oni febst_oni marst_oni aprst_oni mayst_oni junst_oni julst_oni augst_oni sepst_oni octst_oni novst_oni decst_oni]

    global ts_oni_st = zeros(tlength)
    for i in 1:12:tlength
      ts_oni_st[i]     = st_oni[1]
      ts_oni_st[i+1]   = st_oni[2]
      ts_oni_st[i+2]   = st_oni[3]
      ts_oni_st[i+3]   = st_oni[4]
      ts_oni_st[i+4]   = st_oni[5]
      ts_oni_st[i+5]   = st_oni[6]
      ts_oni_st[i+6]   = st_oni[7]
      ts_oni_st[i+7]   = st_oni[8]
      ts_oni_st[i+8]   = st_oni[9]
      ts_oni_st[i+9]   = st_oni[10]
      ts_oni_st[i+10]  = st_oni[11]
      ts_oni_st[i+11]  = st_oni[12]
    end
    # ts_oni_st is the standard deviation of ONI: sigma(oni)
#-------------------------------------------------------------------------- 


    ss1  = zeros(12)
    jan = [mean([tropmn_ts[i] for i in 1:12:tlength])]
    feb = [mean([tropmn_ts[i] for i in 2:12:tlength])]
    mar = [mean([tropmn_ts[i] for i in 3:12:tlength])]
    apr = [mean([tropmn_ts[i] for i in 4:12:tlength])]
    may = [mean([tropmn_ts[i] for i in 5:12:tlength])]
    jun = [mean([tropmn_ts[i] for i in 6:12:tlength])]
    jul = [mean([tropmn_ts[i] for i in 7:12:tlength])]
    aug = [mean([tropmn_ts[i] for i in 8:12:tlength])]
    sep = [mean([tropmn_ts[i] for i in 9:12:tlength])]
    oct = [mean([tropmn_ts[i] for i in 10:12:tlength])]
    nov = [mean([tropmn_ts[i] for i in 11:12:tlength])]
    dec = [mean([tropmn_ts[i] for i in 12:12:tlength])]
   
    ss1 = [jan feb mar apr may jun jul aug sep oct nov dec]

    # compute ts for the tropical mean sst minus seasonal cycle
    global ts_trmn_nsc = zeros(tlength)
    for i in 1:12:tlength
      ts_trmn_nsc[i]    = tropmn_ts[i]    - ss1[1]
      ts_trmn_nsc[i+1]  = tropmn_ts[i+1]  - ss1[2]
      ts_trmn_nsc[i+2]  = tropmn_ts[i+2]  - ss1[3]
      ts_trmn_nsc[i+3]  = tropmn_ts[i+3]  - ss1[4]
      ts_trmn_nsc[i+4]  = tropmn_ts[i+4]  - ss1[5]
      ts_trmn_nsc[i+5]  = tropmn_ts[i+5]  - ss1[6]
      ts_trmn_nsc[i+6]  = tropmn_ts[i+6]  - ss1[7]
      ts_trmn_nsc[i+7]  = tropmn_ts[i+7]  - ss1[8]
      ts_trmn_nsc[i+8]  = tropmn_ts[i+8]  - ss1[9]
      ts_trmn_nsc[i+9]  = tropmn_ts[i+9]  - ss1[10]
      ts_trmn_nsc[i+10] = tropmn_ts[i+10] - ss1[11]
      ts_trmn_nsc[i+11] = tropmn_ts[i+11] - ss1[12]
    end

    # smooth the time series with a running mean
    global ts_nino34_rmn = zeros(tlength)
    global ts_tropmn_rmn = zeros(tlength)
    global ts_rmn        = zeros(tlength)
    global ts_rmn2       = zeros(tlength)
    global oniMtpm       = zeros(tlength)
    istart= 2
    jend  = tlength
    for i in istart:jend-1
      ts_nino34_rmn[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3;
      ts_tropmn_rmn[i] = (ts_trmn_nsc[i+1]+ts_trmn_nsc[i]+ts_trmn_nsc[i-1])/3;
    end
    ts_nino34_rmn[1]=ts_nino34_rmn[2]
    ts_nino34_rmn[jend]=ts_nino34_rmn[jend-1]
    ts_tropmn_rmn[1]=ts_tropmn_rmn[2]
    ts_tropmn_rmn[jend]=ts_tropmn_rmn[jend-1]

    #--------------------------------------------------
    # compute the unscaled RONI
    oniMtpm = ts_oni_cmip .- ts_tropmn_rmn
    #--------------------------------------------------

    #--------------------------------------------------
    # compute the standard deviation of oni - tropmn
    st_oniMtm    = zeros(12)
    janst_oniMtm = std(oniMtpm[i] for i in 1:12:tlength)
    febst_oniMtm = std(oniMtpm[i] for i in 2:12:tlength)
    marst_oniMtm = std(oniMtpm[i] for i in 3:12:tlength)
    aprst_oniMtm = std(oniMtpm[i] for i in 4:12:tlength)
    mayst_oniMtm = std(oniMtpm[i] for i in 5:12:tlength)
    junst_oniMtm = std(oniMtpm[i] for i in 6:12:tlength)
    julst_oniMtm = std(oniMtpm[i] for i in 7:12:tlength)
    augst_oniMtm = std(oniMtpm[i] for i in 8:12:tlength)
    sepst_oniMtm = std(oniMtpm[i] for i in 9:12:tlength)
    octst_oniMtm = std(oniMtpm[i] for i in 10:12:tlength)
    novst_oniMtm = std(oniMtpm[i] for i in 11:12:tlength)
    decst_oniMtm = std(oniMtpm[i] for i in 12:12:tlength)
    st_oniMtm    = [janst_oniMtm febst_oniMtm marst_oniMtm aprst_oniMtm mayst_oniMtm junst_oniMtm julst_oniMtm augst_oniMtm sepst_oniMtm octst_oniMtm novst_oniMtm decst_oniMtm]

    global tsoniMtm_st = zeros(tlength)
    for i in 1:12:tlength
      tsoniMtm_st[i]     = st_oniMtm[1]
      tsoniMtm_st[i+1]   = st_oniMtm[2]
      tsoniMtm_st[i+2]   = st_oniMtm[3]
      tsoniMtm_st[i+3]   = st_oniMtm[4]
      tsoniMtm_st[i+4]   = st_oniMtm[5]
      tsoniMtm_st[i+5]   = st_oniMtm[6]
      tsoniMtm_st[i+6]   = st_oniMtm[7]
      tsoniMtm_st[i+7]   = st_oniMtm[8]
      tsoniMtm_st[i+8]   = st_oniMtm[9]
      tsoniMtm_st[i+9]   = st_oniMtm[10]
      tsoniMtm_st[i+10]  = st_oniMtm[11]
      tsoniMtm_st[i+11]  = st_oniMtm[12]
    end
    # tsoniMtm_st is the standard deviation of ONI-tropmn: sigma(oni-tmn)
    #--------------------------------------------------

    # compute the scaling factor
    sig_oni   =  std(ts_nino34_rmn)
    sig_diff  =  std(ts_nino34_rmn-ts_tropmn_rmn)
    sig_scale =  sig_oni/sig_diff

    sigma     =  ts_oni_st./tsoniMtm_st

    #print(ts_nino34_rmn)
    #println("fury and hatred")
    # compute the relative oni index (RONI):
    #ts_rmn = sig_scale.*(ts_nino34_rmn-ts_tropmn_rmn);
    ts_rmn = sigma.*oniMtpm;
    ts_rmn2 = sig_scale.*(ts_nino34_rmn);
    return ts_rmn;
end

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
