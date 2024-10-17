import time as tt

import collections as cx

import numpy  as np


def lenstat_engine(
    indf,
    oligolimit,
    liner):
    '''
    Compute the length statistics of the
    elements and the resulting oligos.

    :: indf
       type - pd.DataFrame
       desc - a pandas DataFrame storing
              annotated oligopool variants and their parts
    :: oligolimit
       type - integer
       desc - maximum oligo length allowed in the oligopool,
              must be 4 or greater
    '''

    # Book-keeping
    t0       = tt.time()        # Start Timer
    intstats = cx.OrderedDict() # Stats Storage
    fraglens = None             # Running Column  Lengths
    minoligolength = 0          # Running Minimum Length
    maxoligolength = 0          # Running Maximum Length
    minspaceavail = None        # Minimum Free Space Available
    maxspaceavail = None        # Maximum Free Space Available

    # Compute Columnwise Contribution
    for idx,col in enumerate(indf.columns):

        # Extract Length Contribution from Current Column
        collens = np.array([len(seq.replace('-', '')) for seq in indf[col]])
        minelementlen = np.min(collens)
        maxelementlen = np.max(collens)

        if fraglens is None:
            fraglens = collens
        else:
            fraglens += collens

        minoligolength = np.min(fraglens)
        maxoligolength = np.max(fraglens)

        # Update Stats
        intstats[idx] = [
            col,
            minelementlen,
            maxelementlen,
            minoligolength,
            maxoligolength,
            ('Yes', 'No')[maxoligolength <= oligolimit]]

        # Show Update
        if minelementlen == maxelementlen:
            liner.send(
                ' Element {}: Occupies {:,} Base Pair(s)'.format(
                    col,
                    minelementlen))
        else:
            liner.send(
                ' Element {}: Occupies {:,} to {:,} Base Pair(s)'.format(
                    col,
                    minelementlen,
                    maxelementlen))

    minspaceavail = oligolimit - maxoligolength
    maxspaceavail = oligolimit - minoligolength

    # Show Time Elapsed
    liner.send(
        '|* Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return (intstats,
        minspaceavail,
        maxspaceavail)