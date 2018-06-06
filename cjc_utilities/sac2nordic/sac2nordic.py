#!/usr/bin/python

##########################SCRIPT INFO##########################################

# Function to read in picked SAC files and output them to seisan
# formatted pick files (nordic format) and miniseed waveforms

# Requires a filename argument of a list file with the full path of
# the input SAC files for one event - will output merged mseed waveform
# and nordic S-file to the local directory

# Calum Chamberlain
# v1.0 - adding feautures to basic system - working
# 09/08/2014

###############################################################################
import sys, os, datetime


def sac2nor(filename,netwrk,user):
    # Import modules
    from obspy import read as obsread
    from obspy import Stream, Trace

    # Open listfile
    listfile=open(filename,'r')

    # Need to open nordic file with the info of the event in the filename...
    # To get info, read in first file in filelist
    contents=listfile.read()
    first_file=contents.split('\n', 1)[0]
    listfile.seek(0)
    # Fi    nd out how many files are to be merged
    i=0
    for line in listfile:
        i+=1
    # Read first file to get necessary info
    st=obsread(first_file)
    starttime=st[0].stats.starttime
    # Nordic file called: dd-hhmm-ssT.SyyyyMM where T is type (L,R,D) and MM is month
    nordd   =str(starttime.day).zfill(2)
    norhh   =str(starttime.hour).zfill(2)
    normm   =str(starttime.minute).zfill(2)
    norss   =str(starttime.second).zfill(2)
    noryyyy =str(starttime.year).zfill(4)
    norMM   =str(starttime.month).zfill(2)
    norT    ='L' # Set type to local by default
    norname=nordd+'-'+norhh+normm+'-'+norss+norT+'.S'+noryyyy+norMM
    norfile=open(norname,'w')
    # Waveform filename is yyyy-MM-dd-hhmm-ss.netwrk_len
    wavname=noryyyy+'-'+norMM+'-'+nordd+'-'+norhh+normm+'-'+norss+'.'+netwrk+'_'\
            +str(i).zfill(3)

    # Write nordic header info
    norfile.write(' '+noryyyy+' '+norMM+nordd+' '+norhh+normm+' '+norss+\
                  '.0 LM                      TST                               1\n')
    norfile.write(' '+wavname+'                                                  6\n')
    timenow=datetime.datetime.now()
    ID=noryyyy+norMM+nordd+norhh+normm+norss
    timenow=str(timenow.year)[2:4]+'-'+str(timenow.month).zfill(2)+'-'\
            +str(timenow.day).zfill(2)+' '+str(timenow.hour).zfill(2)+':'\
            +str(timenow.minute).zfill(2)
    norfile.write(' ACTION:NEW '+timenow+' OP:'+user+' STATUS:               ID:'+\
                  ID+'     I\n')
    norfile.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO SNR'+\
                  ' AR TRES W  DIS CAZ7\n')

    # Read in list file line by line
    listfile.seek(0)
    i=0
    for line in listfile:
        sacfile=line.split('\n', 1)[0]
        st=obsread(sacfile, dtype='int32')
        # print sacfile
        tr=st[0]
        # Get station and timing info
        # Longitude
        stalong =tr.stats.sac.stlo
        # Latitude
        stalat  =tr.stats.sac.stla
        # Elevation
        stael   =tr.stats.sac.stel
        # Station name
        sta     =tr.stats.station
        sta=sta.ljust(5)
        # Channel
        chan    =tr.stats.channel
        if len(chan)==3:
            chan=chan[0]+chan[2]
        chan=chan.ljust(3)
        # Waveform start time
        wavstrt =tr.stats.starttime
        # Waveform end time
        wavendt =tr.stats.endtime
        # Sampling frequency
        wavfreq =tr.stats.sampling_rate

        # Get location info
        # Longitude
        evlong  =tr.stats.sac.evlo
        # Latitude
        evlat   =tr.stats.sac.evla
        # Depth
        evdep   =tr.stats.sac.evdp
        # Magnitude
        evmag   =tr.stats.sac.mag

    ##### Get pick info#####
        pick_time1=tr.stats.sac.a # Look in first probable header location
        pick_time2=tr.stats.sac.t0 # Look in second probable header location
        if pick_time1 != -12345.0:
            # Format the times
            p_time=wavstrt+pick_time1
            p_hh=str(p_time.hour).zfill(2)
            p_mm=str(p_time.minute).zfill(2)
            p_ss=str(p_time.second).zfill(2)+'.'+str(p_time.microsecond)[0:2]
            p_type=tr.stats.sac.ka[0:2]
            p_type=p_type.ljust(5)
            if p_type=='P    ':
                p_type='IP   '
            elif p_type=='S    ':
                p_type='IS   '
            p_polarity=tr.stats.sac.ka[2:3]
            p_weight=tr.stats.sac.ka[3:4]
            if p_weight=='u':
                p_weight='2' # Arbitrarily set unkown quality picks to 2
            print p_type
        elif pick_time2 != -12345.0:
            p_time=wavstrt+pick_time2
            p_hh=str(p_time.hour).zfill(2)
            p_mm=str(p_time.minute).zfill(2)
            p_ss=str(p_time.second).zfill(2)+'.'+str(p_time.microsecond)[0:2]
            p_type=tr.stats.sac.kt0[0:2]
            p_type=p_type.ljust(5)
            if p_type=='P    ':
                p_type='IP   '
            elif p_type=='S    ':
                p_type='IS   '
            p_polarity=tr.stats.sac.kt0[2:3]
            p_weight=tr.stats.sac.kt0[3:4]
            if p_weight=='u':
                p_weight='2' # Arbitrarily set unkown quality picks to 2
            print p_type
        else:
            p_type=[]
        # Play with data to put out the correct waveform data
        wavedata=tr.data
        wavestats={'network':   tr.stats.network, 'station':    sta, 'location': '',
                   'channel':   tr.stats.channel, 'npts':   tr.stats.npts,'sampling_rate':  wavfreq,
                   'mseed': {'dataquality': 'D'}}
        wavestats['starttime']=wavstrt
        stmseed=Stream([Trace(data=wavedata, header=wavestats)])
        # Add waveform to stream to be written out
        if i != 0:
            stout+=stmseed
        else:
            stout=stmseed
        # Write the pick info to the nordic file
        if p_type:
            print(' '+sta+chan+p_type+p_weight+'   '+p_hh+p_mm+' '+p_ss+'                                                    ')
            norfile.write(' '+sta+chan+p_type+p_weight+'   '+p_hh+p_mm+' '+p_ss+'                                                    \n')
        i+=1
    listfile.close()
    # Close the nordic file
    norfile.close()
    # Write the waveform file
    # print stout[0].stats
    stout.write(wavname,format='MSEED', encoding=11, reclen=256)
    print 'Written wavefrom file: '+wavname+' and nordic file: '+norname

if __name__ == '__main__':
    # Read arguments
    if len(sys.argv) != 4:
        print 'Requires three arguments: sac2nordic <filelist> <network> <user>'
        print 'User should be your 4 character name used in seisan'
        sys.exit()
    else:
        filename=str(sys.argv[1])
        netwrk=str(sys.argv[2])
        user=str(sys.argv[3])[0:4]

    # Make netwrk name the right length (5 characters)
    if len(netwrk)<5:
        netwrk=netwrk.ljust(5,'_')
    elif len(netwrk)>5:
        print 'Network name must be 5 characters or less'
        sys.exit()

    # Check filename exists
    if os.path.isfile(filename) == False:
        print 'No list file found at ',filename
        sys.exit()

    sac2nor(filename,netwrk,user)
