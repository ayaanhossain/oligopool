import time as tt

import collections as cx
import itertools   as ix
from typing import Iterable

import utils as ut


# Parser and Setup Functions

def get_parsed_barcodes(
    bardf,
    barcodes,
    liner):
    '''
    Check barcode feasibility.
    Internal use only.

    :: bardf
       type - pd.DataFrame
       desc - pandas DataFrame storing
              barcode and constant info
    :: barcodes
       type - pd.Series
       desc - ordered list of barcodes
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Parse Barcode Feasibility
    liner.send(' Parsing Barcodes ...')

    # Build ID Dictionary
    IDdict = dict(enumerate(bardf.index))

    # Parse Barcodes
    barcodedict    = {}
    minbarcodelen  = float('+inf')
    maxbarcodelen  = float('-inf')
    barcodelenuniq = True
    barcodecount   = 0

    # Duplicate Bookkeeping
    seen = cx.defaultdict(
        lambda: list())
    duplicates = set()

    # Extract Barcodes and Length Info
    for idx,barcode in enumerate(barcodes):
        barcodedict[idx] = barcode
        minbarcodelen = min(minbarcodelen, len(barcode))
        maxbarcodelen = max(maxbarcodelen, len(barcode))
        seen[barcode].append(idx)

    # Collect Duplicates
    for idxs in seen.values():
        if len(idxs) > 1:
            duplicates.add(tuple(map(
                lambda x: IDdict[x], idxs)))

    # Finalize Set and Length Uniqueness
    barcodelenuniq = minbarcodelen == maxbarcodelen
    barcodelen     = maxbarcodelen
    barcodecount   = len(seen)
    barcodesuniq   = not duplicates

    # Show Updates
    liner.send(
        ' Barcode Extracted: {:,} Unique Barcode(s)\n'.format(
            barcodecount))

    barcodemsg = ['No [INFEASIBLE] (Duplicate Barcodes Detected)', 'Yes'][barcodesuniq]
    liner.send(
        ' Barcode Unique   : {}\n'.format(
            barcodemsg))

    if barcodelenuniq:
        liner.send(
            ' Barcode Length   : {:,} Base Pair(s)\n'.format(
                barcodelen))
    else:
        liner.send(
            ' Barcode Length   : {:,} to {:,} bp [INFEASIBLE] (Variable Length Barcodes)\n'.format(
                minbarcodelen,
                maxbarcodelen))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Compute parsestatus
    parsestatus = all([
        barcodelenuniq,
        barcodesuniq])

    # Show Final Update
    if parsestatus:
        liner.send(
            ' Verdict: Indexing Possibly Feasible\n')
    else:
        liner.send(
            ' Verdict: Indexing Infeasible due to Input Barcode(s) and/or Constant(s)\n')

    return (parsestatus,
        IDdict,
        barcodedict,
        duplicates,
        barcodecount,
        barcodesuniq,
        barcodelen,
        barcodelenuniq,
        minbarcodelen,
        maxbarcodelen)

def get_trimmed_constant(
    constant,
    targetlen,
    constanttype):
    '''
    Return trimmed constant.
    Internal use only.

    :: constant
       type - string
       desc - constant to trim
    :: targetlen
       type - integer
       desc - maximum constant length
              after trimming
    :: constanttype
       type - integer
       desc - 0 = prefix constant
              1 = suffix constant
    '''

    # Extract Suffix from Prefix
    if constanttype == 0:
        return constant[-targetlen:]
    # Extract Prefix from Suffix
    else:
        return constant[:targetlen]

def get_parsed_constants(
    prefixconstant,
    suffixconstant,
    liner):
    '''
    Check barcode constant feasibility.
    Internal use only.

    :: prefixconstant
       type - pd.Series / None
       desc - prefix constant to barcodes
    :: suffixconstant
       type - pd.Series / None
       desc - suffix constant to barcodes
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Parse Anchor Constant Feasibility
    liner.send(' Parsing Constants ...')

    constantscount = 0
    constantsuniq  = None
    longconstants  = None
    prefixuniq     = None
    suffixuniq     = None
    prefixlen      = None
    suffixlen      = None

    # Parse Prefix Constant
    if not prefixconstant is None:
        constantscount += 1
        prefixuniq      = ut.get_unique_count(
            iterable=prefixconstant) == 1

        if constantsuniq is None:
            constantsuniq = True
        constantsuniq = constantsuniq and prefixuniq

        if prefixuniq:
            prefixconstant = get_trimmed_constant(
                constant=prefixconstant[0],
                targetlen=20,
                constanttype=0)
            prefixlen = len(prefixconstant)

            if longconstants is None:
                longconstants = True
            longconstants = longconstants and prefixlen >= 6
        else:
            prefixconstant = None

    # Parse Suffix Constant
    if not suffixconstant is None:
        constantscount += 2
        suffixuniq      = ut.get_unique_count(
            iterable=suffixconstant) == 1

        if constantsuniq is None:
            constantsuniq = True
        constantsuniq = constantsuniq and suffixuniq

        if suffixuniq:
            suffixconstant = get_trimmed_constant(
                constant=suffixconstant[0],
                targetlen=20,
                constanttype=1)
            suffixlen = len(suffixconstant)

            if longconstants is None:
                longconstants = True
            longconstants = longconstants and suffixlen >= 6
        else:
            suffixconstant = None

    # Are the Anchors Feasible?
    constantsextracted = constantscount > 0
    anchormsg = {
        0: 'None Available [INFEASIBLE] (At Least One Constant Required)',
        1: 'Only Prefix',
        2: 'Only Suffix',
        3: 'Both Prefix and Suffix'
    }[constantscount]

    # Show Updates
    liner.send(
        ' Constants Extracted: {}\n'.format(
            anchormsg))

    if not prefixuniq is None:
        uniqmsg = ['No [INFEASIBLE] (Multiple Prefixes Found)', 'Yes'][prefixuniq]
        liner.send(
            '    Prefix Unique   : {}\n'.format(
                uniqmsg))
    if not prefixlen is None:
        uniqmsg = [' [INFEASIBLE] (Prefixes Shorter than 6 bp)', ''][prefixlen >= 6]
        liner.send(
            '    Prefix Length   : {:,} bp{}\n'.format(
                prefixlen,
                uniqmsg))

    if not suffixuniq is None:
        uniqmsg = ['No [INFEASIBLE] (Multiple Suffixes Found)', 'Yes'][suffixuniq]
        liner.send(
            '    Suffix Unique   : {}\n'.format(
                uniqmsg))
    if not suffixlen is None:
        uniqmsg = [' [INFEASIBLE] (Suffixes Shorter than 6 bp)', ''][suffixlen >= 6]
        liner.send(
            '    Suffix Length   : {:,} bp{}\n'.format(
                suffixlen,
                uniqmsg))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Compute parsestatus
    parsestatus = all([
        constantsextracted,
        constantsuniq,
        longconstants])

    # Show Final Update
    if parsestatus:
        liner.send(
            ' Verdict: Indexing Possibly Feasible\n')
    else:
        liner.send(
            ' Verdict: Indexing Infeasible due to Barcode Constant(s)\n')

    # Return Results
    return (parsestatus,
        constantsextracted,
        constantsuniq,
        longconstants,
        prefixconstant,
        suffixconstant,
        prefixuniq,
        suffixuniq,
        prefixlen,
        suffixlen)

def get_extracted_variants(
    variants,
    expcount,
    IDdict,
    warn,
    liner):
    '''
    Return extracted variant sequences.
    Internal use only.

    :: variants
       type - pd.Series
       desc - ordered list of barcode
              associated variants
    :: expcount
       type - integer
       desc - total number of variants
              expected based on the
              barcode count
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Extract Variants
    liner.send(' Parsing Variants ...')

    variantdict  = dict(enumerate(variants))
    variantcount = ut.get_unique_count(
        iterable=variants)
    variantuniq  = variantcount == expcount

    # Compute Duplicates
    duplicates = set()
    if not variantuniq:
        seen = cx.defaultdict(
            lambda: list())
        for idx,variant in variantdict.items():
            seen[variant].append(idx)
        for variant,idxs in seen.items():
            if len(idxs) > 1:
                duplicates.add(tuple(map(
                    lambda x: IDdict[x], idxs)))

    # Show Updates
    liner.send(
        ' Variant Extracted: {:,} Unique Variant(s)\n'.format(
            variantcount))

    variantmsg = ['No [WARNING] (Duplicate Variants Detected)', 'Yes'][variantuniq]
    liner.send(
        ' Variant Unique   : {}\n'.format(
            variantmsg))
    if not variantuniq:
        warn['warncount'] += len(duplicates)
        warn['vars'] = {
            'expectedcount': expcount,
             'variantcount': variantcount,
            'variantunique': variantuniq,
               'duplicates': duplicates}

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return variantdict

# Engine Helper Functions

def infer_kvalue(
    elementlen,
    minimum):
    '''
    Infer k-values for barcode set.
    Internal use only.

    :: elementlen
       type - iterable
       desc - an list of reference
              sequences
    :: minimum
       type - integer
       desc - minimum k-value
    '''

    return max(minimum, (elementlen // 2) - 1)

def infer_tvalue(
    elementlen,
    maximum):
    '''
    Infer t-values for index element.
    Internal use only.

    :: elementlen
       type - iterable
       desc - an list of reference
              sequences
    :: minimum
       type - integer
       desc - maximum t-value
    '''

    return min(maximum,
              (elementlen - 1) // 10)

def split_map(
    maptosplit,
    seqstr=0,
    seqlen=None):
    '''
    Split a map in to key-value
    arrays. Internal use only.

    :: maptosplit
       type - dict
       desc - a dictionary to split in
              to key-value arrays
    :: seqstr
       type - integer
       desc - sequence prefix starting
              index
              (default=0)
    :: seqlen
       type - None / integer
       desc - sequence prefix length
              (default=None)
    '''

    X, Y = [], []
    while maptosplit:
        y,x = maptosplit.popitem()
        X.append(x[seqstr:seqlen])
        Y.append(y)
    return X, Y