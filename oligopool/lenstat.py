import time as tt

import collections as cx
import numpy       as np
import pandas      as pd

from .base import utils    as ut
from .base import valparse as vp


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

def lenstat(
    indata,
    oligolimit,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Length Statistics Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Length Statistics]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field=' Input Data ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field=' Oligo Limit',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)

    # Length Stats Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Launching Length Stats Computation
    liner.send('\n[Step 1: Computing Length Statistics]\n')

    # Compute Length Stats
    (intstats,
    minspaceavail,
    maxspaceavail) = lenstat_engine(
        indf=indf,
        oligolimit=oligolimit,
        liner=liner)

    # Length Statistics Results
    liner.send('\n[Length Statistics]\n')

    # Build Stats Print String
    statsprint = '\n'.join(
        ' ' + line for line in pd.DataFrame.from_dict(
            {k: [u[0]] + [str(v) + ' bp' for v in u[1:-1]] + [u[-1]] \
                for k,u in intstats.items()},
            orient='index',
            columns=pd.MultiIndex.from_arrays(
                [('', '    Min', '    Max', '   Min', '   Max', '   Oligo'),
                 ('', 'Element', 'Element', ' Oligo', ' Oligo', '   Limit'),
                 ('', ' Length', ' Length', 'Length', 'Length', 'Overflow')]
                )).to_string(index=False).split('\n'))

    # Show Stats String
    print('\n{}'.format(
        statsprint))

    # Restructure Length Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'computing-length-statistics',
        'vars'    : {
            'len_stat'     : cx.OrderedDict({v[0]: { # Store Element-wise Length Stats
                'min_element_len':  v[1],
                'max_element_len':  v[2],
                'min_oligo_length': v[3],
                'max_oligo_length': v[4],
                'limit_overflow':   v[5] == 'Yes'} for k,v in intstats.items()}),
            'oligo_limit'    : oligolimit,            # Specified Oligo Limit
            'min_space_avail': minspaceavail,         # Minimum Space Available
            'max_space_avail': maxspaceavail},        # Maximum Space Available
        'warns'   : warns}

    # Show Free Space Available
    if minspaceavail == maxspaceavail:
        liner.send(
            '\n Free Space Available: {:,} Base Pair(s)\n'.format(
                minspaceavail))
    else:
        liner.send(
            '\n Free Space Available: {:,} to {:,} Base Pair(s)\n'.format(
                minspaceavail,
                maxspaceavail))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Close Liner
    liner.close()

    # Return Results
    return stats
