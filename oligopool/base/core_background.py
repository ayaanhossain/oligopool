import time as tt

from . import vectordb as db
from . import utils as ut


# Background Feasibility Check

def is_background_feasible(
    seq,
    bg,
    leftcontext,
    rightcontext):
    '''
    Check if sequence avoids background
    k-mers including junction regions.
    Internal use only.

    :: seq
       type - string
       desc - designed element sequence
    :: bg
       type - vectorDB / list / None
       desc - background k-mer database(s)
    :: leftcontext
       type - string
       desc - left flanking sequence
    :: rightcontext
       type - string
       desc - right flanking sequence
    '''

    # No Background?
    if bg is None:
        return True

    # Empty Element?
    if len(seq) == 0:
        return True

    # Normalize to list
    backgrounds = bg if isinstance(bg, list) else [bg]

    # Must pass ALL backgrounds
    for background in backgrounds:

        # Skip None entries
        if background is None:
            continue

        # Trim Context to Relevant Region
        # Only need k-1 bases from each side for junction k-mers
        leftctx  = leftcontext[-(background.K-1):]  if leftcontext  and len(leftcontext)  >= background.K-1 else (leftcontext  or '')
        rightctx = rightcontext[:background.K-1]    if rightcontext and len(rightcontext) >= background.K-1 else (rightcontext or '')
        fullseq  = leftctx + seq + rightctx

        # Compute Check Region
        # First k-mer: starts where it first touches element
        # Last k-mer: starts at last element base
        elemstart = len(leftctx)
        elemend   = len(leftctx) + len(seq) - 1

        # First k-mer touching element starts at max(0, elemstart - K + 1)
        # Last k-mer touching element starts at elemend
        start = max(0, elemstart - background.K + 1)
        end   = elemend + 1  # exclusive

        # Check All k-mers in Region
        for i in range(start, end):
            kmer = fullseq[i:i+background.K]
            if len(kmer) == background.K and kmer in background:
                return False

    # No Conflicts
    return True


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
    stats['vars']['filled_kmer_count'] = fill_count
    stats['vars']['remaining_kmer_count'] = left_count

    # If Successful Update and Close DB
    if stats['status']:
        vDB.DB['LEN'] = fill_count
        vDB.close()

    # Otherwise Drop DB
    else:
        vDB.drop()

    # Return Stats
    return stats
