import time as tt

import collections as cx

from . import utils as ut
from . import scry  as sy


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
        ' Barcode Extracted: {:,} Unique Sequences(s)\n'.format(
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
    attachetype,
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
    :: attachetype
       type - integer
       desc - constant attachement type
              0 = barcode
              1 = associate
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
                constant=prefixconstant.iloc[0],
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
                constant=suffixconstant.iloc[0],
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
        attache = ['Barcode', 'Associate'][attachetype]
        liner.send(
            ' Verdict: Indexing Infeasible due to {} Constant(s)\n'.format(
                attache))

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

def get_extracted_associates(
    associates,
    expcount,
    IDdict,
    warn,
    liner):
    '''
    Return extracted associate sequences.
    Internal use only.

    :: associates
       type - pd.Series
       desc - ordered list of barcode
              associated associates
    :: expcount
       type - integer
       desc - total number of associates
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

    # Extract Associates
    liner.send(' Parsing Associates ...')

    associatedict  = dict(enumerate(associates))
    associatecount = ut.get_unique_count(
        iterable=associates)
    associateuniq  = associatecount == expcount

    # Compute Duplicates
    duplicates = set()
    if not associateuniq:
        seen = cx.defaultdict(
            lambda: list())
        for idx,associate in associatedict.items():
            seen[associate].append(idx)
        for associate,idxs in seen.items():
            if len(idxs) > 1:
                duplicates.add(tuple(map(
                    lambda x: IDdict[x], idxs)))

    # Show Updates
    liner.send(
        ' Associate Extracted: {:,} Unique Sequence(s)\n'.format(
            associatecount))

    associatemsg = ['No [WARNING] (Duplicate Associates Detected)', 'Yes'][associateuniq]
    liner.send(
        ' Associate Unique   : {}\n'.format(
            associatemsg))
    if not associateuniq:
        warn['warncount'] += len(duplicates)
        warn['vars'] = {
              'expectedcount': expcount,
             'associatecount': associatecount,
            'associateunique': associateuniq,
                 'duplicates': duplicates}

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return associatedict

# Engine Helper Functions

def infer_tvalue(
    elementlen,
    maximum):
    '''
    Infer t-values for index element.
    Internal use only.

    :: elementlen
       type - integer
       desc - length of element
    :: maximum
       type - integer
       desc - maximum t-value
    '''

    return min(maximum, ut.get_tvalue(
        elementlen=elementlen))

def infer_kvalue(
    elementlen,
    tvalue,
    minimum):
    '''
    Infer k-values for barcode set.
    Internal use only.

    :: elementlen
       type - integer
       desc - length of element
    :: t-value
       type - estimated t-value
              for barcode set
    :: minimum
       type - integer
       desc - minimum k-value
    '''

    return max(minimum, (elementlen // (tvalue * 2)))

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

def get_associate_tvalues(
    associatedict,
    maximum,
    liner):
    '''
    Return a dictionary of associate
    t-values. Internal use only.

    :: associatedict
       type - dict
       desc - complete associate sequence
              dictionary
    :: minimum
       type - integer
       desc - maximum t-value
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    associatetvaldict = {}
    mintval = None
    maxtval = None
    for ID,associate in associatedict.items():
        tvalue = infer_tvalue(
            elementlen=len(associate),
            maximum=maximum)
        liner.send(
            ' Inferring Associate t-value: Associate {:,} w/ {:,} Mismatch(es)'.format(
                ID,
                tvalue))
        associatetvaldict[ID] = tvalue
        mintval = min(mintval, tvalue) if not mintval is None else tvalue
        maxtval = max(maxtval, tvalue) if not maxtval is None else tvalue
    return associatetvaldict, mintval, maxtval

def index_engine(
    IDdict,
    barcodedict,
    barcodename,
    barcodecount,
    barcodelen,
    barcodeprefix,
    barcodesuffix,
    barcodepregap,
    barcodepostgap,
    associatedict,
    associateprefix,
    associatesuffix,
    indexdir,
    indexqueue,
    liner):
    '''
    Prepare final indexes and models for
    counting, based on parsed objects.
    Internal use only.

    :: IDdict
       type - dict
       desc - barcode/associate ID dictionary
    :: barcodedict
       type - dict
       desc - complete barcode sequence
              dictionary
    :: barcodename
       type - string
       desc - barcode column name
    :: barcodecount
       type - integer
       desc - total number of barcodes
    :: barcodelen
       type - integer
       desc - length of barcodes
    :: barcodeprefix
       type - string / None
       desc - constant barcode prefix
    :: barcodesuffix
       type - string / None
       desc - constant barcode suffix
    :: barcodepregap
       type - integer
       desc - gap between prefix and barcode
    :: barcodepostgap
       type - integer
       desc - gap between barcode and suffix
    :: associatedict
       type - dict
       desc - complete associate sequence
              dictionary
    :: associateprefix
       type - string / None
       desc - constant associate prefix
    :: associatesuffix
       type - string / None
       desc - constant associate suffix
    :: indexdir
       type - string
       desc - path to directory temporarily
              storing indexed structures
              and data models
    :: indexqueue
       type - mp.SimpleQueue
       desc - queue storing indexed object
              paths once saved into indexdir
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Base Index Directory Path
    indexdir += '/'

    # Counting Meta Data Initialization
    maxtval = 2
    metamap = {
        'variantcount'   : barcodecount,
        'barcodename'    : barcodename,
        'barcodelen'     : barcodelen,
        'barcodeprefix'  : barcodeprefix,
        'barcodesuffix'  : barcodesuffix,
        'barcodepregap'  : barcodepregap,
        'barcodepostgap' : barcodepostgap,
        'barcodegapped'  : any((barcodepregap, barcodepostgap)),
        'association'    : len(associatedict) > 0,
        'associateprefix': associateprefix,
        'associatesuffix': associatesuffix
    }

    # Compute Barcode Set t-value
    liner.send(' Inferring Barcode t-value ...')
    barcodetval = infer_tvalue(
        elementlen=barcodelen,
        maximum=2)
    liner.send(
        '   Barcode t-value: {:,} Mismatch(es)\n'.format(
            barcodetval))
    metamap['barcodetval'] = barcodetval

    # Compute Barcode Set k-value
    liner.send(' Inferring Barcode k-value ...')
    barcodekval = infer_kvalue(
        elementlen=barcodelen,
        tvalue=barcodetval,
        minimum=3)
    liner.send(
        '   Barcode k-value: {:,} Base Pair(s)\n'.format(
            barcodekval))
    metamap['barcodekval'] = barcodekval

    # Compute Associate Set t-value
    if associatedict:
        liner.send(' Inferring Associate t-value ...')
        (associatetvaldict,
        mintval,
        maxtval) = get_associate_tvalues(
            associatedict=associatedict,
            maximum=4,
            liner=liner)
        if mintval == maxtval:
            liner.send(
                ' Associate t-value: {:,} Mismatch(es)\n'.format(
                    mintval))
        else:
            liner.send(
                ' Associate t-value: {:,} to {:,} Mismatch(es)\n'.format(
                    mintval,
                    maxtval))
        metamap['associatetvalmin'] = mintval
        metamap['associatetvalmax'] = maxtval
    else:
        metamap['associatetvalmin'] = None
        metamap['associatetvalmax'] = None

    # Compute Barcode Prefix t-value
    if not barcodeprefix is None:
        bpxtval = infer_tvalue(
            elementlen=len(barcodeprefix),
            maximum=2)
        metamap['bpxtval'] = bpxtval
    else:
        metamap['bpxtval'] = None

    # Compute Barcode Suffix t-value
    if not barcodesuffix is None:
        bsxtval = infer_tvalue(
            elementlen=len(barcodesuffix),
            maximum=2)
        metamap['bsxtval'] = bsxtval
    else:
        metamap['bsxtval'] = None

    # Compute Trimming Policy
    metamap['trimpolicy'] = None
    if not barcodeprefix is None and \
       not barcodesuffix is None:
        metamap['trimpolicy'] = 1.5
    elif not barcodeprefix is None:
        metamap['trimpolicy'] = 1
    elif not barcodesuffix is None:
        metamap['trimpolicy'] = 2

    # Compute Read Anchor Motif
    if not barcodeprefix is None:
        metamap['anchor']     = barcodeprefix
        metamap['revanchor']  = ut.get_revcomp(
            seq=barcodeprefix)
        metamap['anchortval'] = bpxtval
    else:
        metamap['anchor']     = barcodesuffix
        metamap['revanchor']  = ut.get_revcomp(
            seq=barcodesuffix)
        metamap['anchortval'] = bsxtval

    # Compute Associate Prefix t-value
    if not associateprefix is None:
        apxtval = infer_tvalue(
            elementlen=len(associateprefix),
            maximum=2)
        metamap['apxtval'] = apxtval
    else:
        metamap['apxtval'] = None

    # Compute Associate Suffix t-value
    if not associatesuffix is None:
        asxtval = infer_tvalue(
            elementlen=len(associatesuffix),
            maximum=2)
        metamap['asxtval'] = asxtval
    else:
        metamap['asxtval'] = None

    # Save and Delete ID Dict
    liner.send(' Writing ID Map ...')
    opath = indexdir+'ID.map'
    ut.savedict(
        dobj=IDdict,
        filepath=opath)
    indexqueue.put(opath)
    del IDdict
    liner.send(' Writing        ID Map  : Done\n')

    # Save and Delete metamap
    liner.send(' Writing Meta Map ...')
    opath = indexdir+'meta.map'
    ut.savedict(
        dobj=metamap,
        filepath=opath)
    indexqueue.put(opath)
    del metamap
    liner.send(' Writing      Meta Map  : Done\n')

    # Train Barcode Model
    liner.send(' Training Barcode Model ...')
    X,Y = split_map(
        maptosplit=barcodedict)
    t0  = tt.time()
    barcodemodel = sy.Scry().train(
        X=X,
        Y=Y,
        n=barcodelen,
        k=barcodekval,
        t=barcodetval,
        liner=liner)

    # Save Barcode Model
    opath = indexdir+'barcode.model'
    ut.savemodel(
        mobj=barcodemodel,
        filepath=opath)
    indexqueue.put(opath)
    del barcodedict
    del barcodemodel
    liner.send(' Writing   Barcode Model: Done\n')

    # Update Associate Dictionary
    if associatedict:
        liner.send(' Updating Associate Map ...')
        for k in associatedict:
            associatedict[k] = (
                associatedict[k],
                associatetvaldict.pop(k))

    # Save and Delete associatedict
    if associatedict:
        liner.send(' Writing Associate Map ...')
        opath = indexdir+'associate.map'
        ut.savedict(
            dobj=associatedict,
            filepath=opath)
        indexqueue.put(opath)
        del associatedict
        liner.send(' Writing Associate Map  : Done\n')

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Indexing Completed
    indexqueue.put(None)