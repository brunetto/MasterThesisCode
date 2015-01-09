import numpy as np
import sys

class head:
    def __init__(self, fname):
        import types
        bo_mark='>'
        # start by reading in header:    
        if type(fname) is types.StringType:
            f=open(fname, 'rb')
        elif type(fname) is not types.FileType:
            raise TypeError('argument must either be an open file or ' + 
                            'a string containing a file name')
        else:
            f=fname

        self.pad = np.fromfile(f, count=1, dtype=bo_mark+'i4')
        # npart is an array containing the number of particles in the
        # file divided by type (gas, ...)
        self.npart = np.fromfile(f, count=6, dtype=bo_mark+'i4')
        self.massarr = np.fromfile(f, count=6, dtype=bo_mark+'f8')
        self.aaa = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.redshift = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.flag_sfr = np.fromfile(f, count=1, dtype=bo_mark+'i4')
        self.flag_feedback = np.fromfile(f, count=1, dtype=bo_mark+'i4')
        self.nall = np.fromfile(f, count=6, dtype=bo_mark+'i4')
        self.cooling_flag = np.fromfile(f, count=1, dtype=bo_mark+'i4')
        self.numfiles = np.fromfile(f, count=1, dtype=bo_mark+'i4')
        self.boxsize = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.Omega = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.OmegaL0 = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.Hubblepar = np.fromfile(f, count=1, dtype=bo_mark+'f8')
        self.version = np.fromfile(f, count=1, dtype=bo_mark+'a96')
        self.pad2=np.fromfile(f, count=1, dtype=bo_mark+'i4')
        if type(fname) is types.StringType: f.close()

def read_gif2_file(fname):
    bo_mark = '>'
    f=open(fname, 'rb')
    # start by reading in header:    
    ghead=head(f)
    npt=ghead.npart.sum()

    f.seek(4, 1)
    pos=np.fromfile(f, count=npt*3, dtype=bo_mark + 'f4').reshape((npt, 3))
    return pos
