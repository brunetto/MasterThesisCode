import numpy as np
import sys

class head:
    """object version of 'header' function.  Holds info about # of
    particles, mass of each particle type, time and redshift, and any
    flags in the header
    
    contents:
    self.npart
    self.mass
    self.time
    self.flags
    self.npartall
    """
    def __init__(self, fname, byteswap=False):
	"""read in the header.
        fname  -- either astring containing the file name or an open file object
        byteswap=False  -- if True, swap the endian-ness on read-in (on little
        endian machines, will assume file is big endian and vice versa)"""
        if byteswap:
            import sys
            # check system byte-ordering.
            byteord=sys.byteorder
            if byteord == 'little':
                bo_mark='<'
            else:
                bo_mark='>'
        else:
            # use native byte-ordering
            bo_mark='='
        import types

        # start by reading in header:    
        if type(fname) is types.StringType:
            f=open(fname, 'rb')
        elif type(fname) is not types.FileType:
            raise TypeError('argument must either be an open file or ' + 
                            'a string containing a file name')
        else:
            f=fname
        pad=np.fromfile(f, count=1, dtype=bo_mark+'i4')
        npart=np.fromfile(f, count=6, dtype=bo_mark+'i4')
        massarr=np.fromfile(f, count=6, dtype=bo_mark+'f8')
        time_reds=np.fromfile(f, count=2, dtype=bo_mark+'f8')
        flag=np.fromfile(f, count=2, dtype=bo_mark+'i4')
        npartall=np.fromfile(f, count=6, dtype=bo_mark+'i4')
        # moreFlags[1] should hold number of files over which snapshot is split
        moreFlags=np.fromfile(f, count=2, dtype=bo_mark+'i4')
        cosmoParams=np.fromfile(f, count=4, dtype=bo_mark+'f8')
        # header is 256 bytes, rest is empty space: 
        empty=np.fromfile(f, count=24, dtype=bo_mark+'i4')
        # done reading header; read f77 record for end of header
        pad=np.fromfile(f, count=1, dtype=bo_mark+'i4')
        if type(fname) is types.StringType: f.close()

	self.npart=npart
	self.mass=massarr
	self.time=time_reds
	self.flags=flag
        self.npartall=npartall
        self.numFiles=moreFlags[1]
        self.cosmoParams=cosmoParams



def read_snapshot(fname, byteswap=False, longIDs=True):
    """
    get data from a GADGET snapshot.  

    INPUTS:
    fname  -- name of GADGET snapshot file
    
    OPTIONAL INPUTS:
        
    byteswap=False  -- if True, swap the endian-ness on read-in (on little
        endian machines, will assume file is big endian and vice versa)
    longIDs=False -- if True, assume ids 8-byte ints, not 4-byte ints
    """
    try: 
        [xx == 'subidorder' for xx in fname.split('_')].index(True)
        reshuffled=True
    except ValueError:
        reshuffled=False

    if longIDs:
        myIDtype='i8'
        myIDlen=8
    else: 
        myIDtype='i4'
        myIDlen=4

    if byteswap:
        import sys
        # check system byte-ordering.
        byteord=sys.byteorder
        if byteord == 'little':
            bo_mark='<'
        else:
            bo_mark='>'
    else:
        # use native byte-ordering
        bo_mark='='

    #if outdat == 'mass':
	#if max(massar) > 0:
         #   raise TypeError('there is no mass block! Exiting')

    f=open(fname, 'rb')
    # start by reading in header:    
    ghead=head(f, byteswap=byteswap)
    npt=ghead.npart.sum() #number of particles???
    #if outdat == 'mass':
	#if ghead.mass.max() > 0:
         #   raise TypeError('there is no mass block! Exiting')

    f.seek(4, 1)
    pos=np.fromfile(f, count=npt*3, dtype=bo_mark + 'f4').reshape((npt, 3)) #positions???
#    f.seek(8, 1)
#    vel=np.fromfile(f, count=npt*3, dtype=bo_mark + 'f4').reshape((npt, 3)) #velocities???
#    f.seek(8, 1)
#    ids=np.fromfile(f, count=npt, dtype=bo_mark + myIDtype) #ids???
    
    snap = {'pos': pos}#, 'vel': vel, 'ids': ids, 'numpart': ghead.npart}

    if not reshuffled:
        f.close()
    else:
#        f.seek(8, 1)
#        hsml=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
#        f.seek(8, 1)
#        dens=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
#        snap['hsml'] = hsml
#        snap['dens'] = dens
        if reshuffled:
            f.seek(8, 1)
            disp=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
            snap['disp'] = disp
        f.close()
    
    return snap
