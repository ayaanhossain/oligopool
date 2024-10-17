import time as tt

from . import vectordb as db
from . import utils as ut


# Engine Objective and Helper Functions

def background_engine(
    background,
    maxreplen,
    outdir,
    stats,
    liner):
    '''
    Extract and populate background.
    Internal use only.

    :: background
       type - list
       desc - list of background sequences
    :: maxreplen
       type - integer
       desc - maximum shared repeat length
    :: outdir
       type - string
       desc - output directory
    :: stats
       type - dict
       desc - background stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Open vectorDB instance
    vDB = db.vectorDB(
        path=outdir,
        maximum_repeat_length=maxreplen)

    # Book-keeping
    t0   = tt.time()
    plen = ut.get_printlen(
        value=len(background))

    # Loop and insert background
    for idx,seq in enumerate(background):

        # Format sequence
        if len(seq) > 20:
            pseq = seq[:20]
            pbuf = '...'
        else:
            pseq = seq
            pbuf = ''

        # Insert sequence
        vDB.add(
            seq=seq,
            rna=False)

        # Show updates
        liner.send(
            ' Sequence {:{},d}: {}{} Inserted'.format(
                idx+1, plen, pseq, pbuf))

    # Final Update
    liner.send(
        ' Sequence {:{},d}: {}{} Inserted\n'.format(
            idx+1, plen, pseq, pbuf))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Populate Stats
    kmer_space = ((4**(maxreplen+1)) // 2)
    fill_count = min(kmer_space, len(vDB))
    left_count = kmer_space - fill_count

    stats['status'] = (left_count * 1.) / kmer_space > .01
    stats['basis']  = 'solved' if stats['status'] else 'infeasible'
    stats['vars']['kmer_space'] = kmer_space
    stats['vars']['fill_count'] = fill_count
    stats['vars']['left_count'] = left_count

    # If Successful Update and Close DB
    if stats['status']:
        vDB.DB.put(
            b'LEN',
            str(fill_count).encode())
        vDB.close()

    # Otherwise Drop DB
    else:
        vDB.drop()

    # Return Stats
    return stats