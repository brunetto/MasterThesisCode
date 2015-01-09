import numpy as np
from numpy import fromfile as np_ff

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

    OUTPUTS
    =======

    Returns outdat(???) or a dictionary.

    """
    outdat = None
    if outdat is not None:
        if outdat not in ['pos', 'vel', 'ids', 'mass', 'hsml', 'dens', 'disp']:
            raise ValueError('outdat must be one of ' +
                             '[pos, vel, ids, mass, hsml, dens, disp]')
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

    if outdat == 'mass':
	if max(massar) > 0:
            raise TypeError('there is no mass block! Exiting')

    f=open(fname, 'rb')
    # start by reading in header:    
    ghead=head(f, byteswap=byteswap)
    npt=ghead.npart.sum()
    if outdat == 'mass':
	if ghead.mass.max() > 0:
            raise TypeError('there is no mass block! Exiting')

    f.seek(4, 1)
    pos=np.fromfile(f, count=npt*3, dtype=bo_mark + 'f4').reshape((npt, 3))
    f.seek(8, 1)
    vel=np.fromfile(f, count=npt*3, dtype=bo_mark + 'f4').reshape((npt, 3))
    f.seek(8, 1)
    ids=np.fromfile(f, count=npt, dtype=bo_mark + myIDtype)
    if not reshuffled:
        f.close()
        if outdat is not None:
            return locals()[outdat]
        else:
            #return pos, vel, ids
            snap  = {'pos': pos, 'vel': vel, 'ids': ids, 'numpart': ghead.npart}
            return snap
    else:
        f.seek(8, 1)
        hsml=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
        f.seek(8, 1)
        dens=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
        if reshuffled:
            f.seek(8, 1)
            disp=np.fromfile(f, count=npt, dtype=bo_mark + 'f4')
        f.close()
        if outdat is not None:
            return locals()[outdat]
        else:
            #return pos, vel, ids, hsml, dens, disp
            snap = {'pos': pos, 'vel': vel, 'ids': ids, 
                           'numpart': ghead.npart, 'dens': dens, 
                           'smooth_lenght': hsml}
            return snap


def read_subtab_file(snapnum, fileNum, longIDs=True, 
                     flag_group_velDisp=False, 
                     only_groups=False, only_subhalos=False, 
                     snapbase="", 
                     return_info=False, quiet=False):
    """
    use for reading in one individual subhalo_tab file

    OUTPUTS
    =======

    Returns a dictionary.

    """
    
    # use dictionaries:
    group={}
    subhalo={}

    if longIDs:
        id_type='int64'
    else:
        id_type='int32'

    # if there are more than 1000 snapshots, need to change 3 -> 4
    dirnum='%03d' % snapnum
#    fname='%s/groups_%s/subhalo_tab_%s.%d' % (snapbase, dirnum, 
#                                              dirnum, fileNum)
    fname = 'subhalo_tab_067.999'

    f=open(fname, 'rb')
    nGr=np_ff(f, dtype='i', count=1)[0]
    tot_nGr=np_ff(f, dtype='i', count=1)[0]
    nIDs=np_ff(f, dtype='i', count=1)[0]
    tot_nIDs=np_ff(f, dtype='int64', count=1)[0]
    nTask=np_ff(f, dtype='i', count=1)[0]
    nSubs=np_ff(f, dtype='i', count=1)[0]
    tot_nSubs=np_ff(f, dtype='i', count=1)[0]
    if return_info:
        f.close()
        info = {'tot_fof_in_file': nGr, 'tot_fof_in_snap': tot_nGr, 
                       'particles_in_fof_in_ID': nIDs, 'particles_in_fof_in_snap': tot_nIDs, 
                       'tot_files_in_snap': nTask, 'sh_in file':nSubs, 'sh_in_snap': tot_nSubs}
        return info

    if nGr > 0:
        group['len'] = np_ff(f, 'i', nGr)
        group['offset'] = np_ff(f, 'u4', nGr)
        group['mass'] = np_ff(f, 'f', nGr)
        group['pos'] = np_ff(f, 'f', nGr*3).reshape(nGr, 3)
        group['m_mean200'] = np_ff(f, 'f', nGr)
        group['r_mean200'] = np_ff(f, 'f', nGr)
        group['m_crit200'] = np_ff(f, 'f', nGr)
        group['r_crit200'] = np_ff(f, 'f', nGr)
        group['m_tophat'] = np_ff(f, 'f', nGr)
        group['r_tophat'] = np_ff(f, 'f', nGr)
        if flag_group_velDisp:
            group['velDisp_mean200'] = np_ff(f, 'f', nGr)
            group['velDisp_crit200'] = np_ff(f, 'f', nGr)
            group['velDisp_tophat'] = np_ff(f, 'f', nGr)
        group['contaminationCount'] = np_ff(f, 'i', nGr)
        group['contaminationMass'] = np_ff(f, 'f', nGr)
        group['n_subs'] = np_ff(f, 'i', nGr)
        group['firstSub'] = np_ff(f, 'i', nGr)

    if not only_groups:
        if nSubs > 0:
            subhalo['len'] = np_ff(f, 'i', nSubs)
            subhalo['offset'] = np_ff(f, 'u4', nSubs)
            subhalo['parent'] = np_ff(f, 'i', nSubs)
            subhalo['mass'] = np_ff(f, 'f', nSubs)
            subhalo['pos'] = np_ff(f, 'f', nSubs*3).reshape((nSubs, 3))
            subhalo['vel'] = np_ff(f, 'f', nSubs*3).reshape((nSubs, 3))
            subhalo['CM'] = np_ff(f, 'f', nSubs*3).reshape((nSubs, 3))
            subhalo['spin'] = np_ff(f, 'f', nSubs*3).reshape((nSubs, 3))
            subhalo['velDisp'] = np_ff(f, 'f', nSubs)
            subhalo['vMax'] = np_ff(f, 'f', nSubs)
            subhalo['rMax'] = np_ff(f, 'f', nSubs)
            subhalo['halfMassRad'] = np_ff(f, 'f', nSubs)
            subhalo['mostBoundID'] = np_ff(f, id_type, nSubs)
            subhalo['GrNr'] = np_ff(f, 'i', nSubs)
    f.close()

    if not quiet:
        print
        print "Total num of groups    =", tot_nGr
        print "Total num of subgroups =", tot_nSubs, '\n'
        print "Number of groups in file %d:     %d" % (fileNum, nGr)
        print "Number of subgroups in file %d:  %d\n" % (fileNum, nSubs)
    if only_groups:
        return group
    elif only_subhalos:
        return subhalo
    else:
        out = {'group': group, 'subhalo': subhalo}
        return out



def read_reshuffled_group(snapnum, groupNum, vel=False, longIDs=True, 
                          outdatType='pos',
                          sfrFlag=False, bhFlag=False, 
                          snapbase='', 
                          snapname='snap_newMillen'):
    """
    get the positions or velocities from a reshuffled snapshot file.
    Note: need to work on subtab files (NOT grouptab files) throughout

    OUTPUT
    ======

    Returns an array?

    """
    #import gadget as g

    # cumulative number of groups read:
    cumNgroups=0
    cumNids=0
    cumNpcls=0
    # tells how many IDs have come before target group, will be used for
    # correcting overflow.
    target_nIDs=0
    # how many particles are in target halo?
    target_np=0

    fileNum=0
    # offsets overflow because of u4:
    maxOffset=2**32
    generalData=read_subtab_file(snapnum, fileNum, longIDs=longIDs, 
                                 snapbase=snapbase, 
                                 quiet=True, return_info=True)
    # number of files:
    nTasks=generalData['particles_in_fof_in_snap']
    # will hold cumulative number of groups for each group_tab file
    lenArr=np.zeros(nTasks, dtype='int32')
    # will hold cumulative number of ids for each snapshot_subidorder file
    idArr=np.zeros(nTasks, dtype='int64')
    for i in range(nTasks):
        groupData=read_subtab_file(snapnum, i, longIDs=longIDs, 
                                     snapbase=snapbase, 
                                     quiet=True, return_info=True)
        nGroups=groupData[0]
        nIDs=groupData[2]
        if cumNgroups <= groupNum:
            nPcls=read_subtab_file(snapnum, i, longIDs=longIDs, 
                                   snapbase=snapbase, quiet=True)[0]['len']
            nPcls_file=nPcls.astype('int64').sum()
            # in this case, target group is in the current file
            if nGroups + cumNgroups > groupNum:
                sum1=nPcls[:groupNum-cumNgroups].astype('int64').sum()
                target_nIDs=cumNpcls + sum1
                target_np=nPcls[groupNum-cumNgroups]
            cumNpcls += nPcls_file
        lenArr[i] = cumNgroups
        idArr[i] = cumNids

        cumNgroups += nGroups
        cumNids += nIDs.astype('int64')
    if cumNgroups != generalData[tot_fof_in_file]:
        raise ValueError('problem with number of groups!')
    if cumNids != generalData[particles_in_fof_in_ID]:
        raise ValueError('problem with number of IDs!')

    
    subFile=lenArr.searchsorted(groupNum, 'right')-1
    # find offset entry of relevant group:
    offset=read_subtab_file(snapnum, subFile, longIDs=longIDs, 
                            snapbase=snapbase, 
                            quiet=True)[0]['offset'][groupNum-lenArr[subFile]]
    
    # correct offset for overflow
    offset += (target_nIDs // maxOffset) * maxOffset
    numSnapFiles=head(snapbase + 'snapdir_%03d/%s_subidorder_%03d.%d' 
                        % (snapnum, snapname, snapnum, 0)).numFiles
    npArr=np.zeros(numSnapFiles, dtype='int64')
    for i in range(numSnapFiles):
        npArr[i]=head(snapbase + 'snapdir_%03d/%s_subidorder_%03d.%d' 
                        % (snapnum, snapname, snapnum, i)).npart[1]

    npArr=npArr.cumsum()
    snapFile1=npArr.searchsorted(offset, 'right')
    snapFile2=npArr.searchsorted(offset + target_np, 'right')


    if outdatType not in ['pos', 'vel', 'ids', 'mass', 'hsml', 'dens', 'disp']:
        raise ValueError('outdat must be one of ' +
                         '[pos, vel, ids, mass, hsml, dens, disp]')
    # handle first file differently from others:
    if offset < npArr[0]:
        startId=offset
    else:
        startId=offset - npArr[snapFile1]
    # startId=offset - npArr[snapFile1]
    # in this case, all IDs are in one file:
#    mysnap='%ssnapdir_%03d/%s_subidorder_%03d.%d' % (snapbase, snapnum, snapname,
#                                                     snapnum, snapFile1)
    mysnap='%s_subidorder_%03d.%d' % (snapbase, snapnum, snapname,
                                                     snapnum, snapFile1)

    print target_np
    if snapFile1 == snapFile2:
        outdat=read_snapshot(mysnap, longIDs=longIDs, 
                               outdat=outdatType)[startId:startId+target_np]

    else:
        # Handles groups split over multiple snapshots; should be working.
        print 'using snapshot files %d-%d' % (snapFile1, snapFile2)
        outdat=np.zeros((target_np,3), dtype='float32')
        firstNum=npArr[snapFile1]-offset
        outdat[:firstNum]=read_snapshot(mysnap, longIDs=longIDs, 
                               outdat=outdatType)[startId:]
        counter=firstNum
        fileToDo=1
        while (fileToDo < snapFile2-snapFile1):
            mysnap='%ssnapdir_%03d/%s_subidorder_%03d.%d' % (snapbase, snapnum, 
                                                             snapname, snapnum, 
                                                             snapFile1+fileToDo)
            numInFile=npArr[snapFile1+fileToDo]-npArr[snapFile1+fileToDo-1]
            outdat[counter:counter+numInFile]=read_snapshot(mysnap, 
                                                              longIDs=longIDs, 
                                                              outdat=outdatType)
            counter += numInFile
            fileToDo += 1
            
        mysnap='%ssnapdir_%03d/%s_subidorder_%03d.%d' % (snapbase, snapnum, 
                                                         snapname,
                                                         snapnum, snapFile2)
        outdat[counter:]=read_snapshot(mysnap, longIDs=longIDs, 
                                         outdat=outdatType)[:target_np-counter]
    
    return outdat


