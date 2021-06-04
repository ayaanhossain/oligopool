import time as tt

import collections     as cx
import itertools       as ix
import multiprocessing as mp
import atexit          as ae

import utils     as ut
import valparse  as vp
import scry      as sy
import coreindex as ci


def index_engine(
    IDdict,
    barcodedict,
    barcodename,
    barcodecount,
    barcodelen,
    barcodeprefix,
    barcodesuffix,
    variantdict,
    variantprefix,
    variantsuffix,
    indexdir,
    indexqueue,
    liner):
    '''
    Prepare final indexes and models for
    counting, based on parsed objects.
    Internal use only.

    :: IDdict
       type - dict
       desc - barcode/variant ID dictionary
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
    :: variantdict
       type - dict
       desc - complete variant sequence
              dictionary
    :: variantprefix
       type - string / None
       desc - constant variant prefix
    :: variantsuffix
       type - string / None
       desc - constant variant suffix
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
        'eventcount'   : barcodecount,
        'barcodename'  : barcodename,
        'barcodelen'   : barcodelen,
        'barcodeprefix': barcodeprefix,
        'barcodesuffix': barcodesuffix,
        'association'  : len(variantdict) > 0,
        'variantprefix': variantprefix,
        'variantsuffix': variantsuffix,
        'IDmap'        : {
            k:IDdict.pop(k) for k in range(len(IDdict))
        }
    }

    # Compute Barcode Set k-value
    liner.send(' Inferring Barcode k-value ...')
    barcodekval = ci.infer_kvalue(
        elementlen=barcodelen,
        minimum=4)
    liner.send(
        ' Barcode k-value: {:,} Base Pair(s)\n'.format(
            barcodekval))
    metamap['barcodekval'] = barcodekval

    # Compute Barcode Set t-value
    liner.send(' Inferring Barcode t-value ...')
    barcodetval = ci.infer_tvalue(
        elementlen=barcodelen,
        maximum=2)
    liner.send(
        ' Barcode t-value: {:,} Mismatch(es)\n'.format(
            barcodetval))
    metamap['barcodetval'] = barcodetval

    # Compute Variant Set t-value
    if variantdict:
        liner.send(' Inferring Barcode t-value ...')
        varianttval = ci.infer_tvalue(
            elementlen=min(map(len, variantdict.values())),
            maximum=2)
        liner.send(
            ' Variant t-value: {:,} Mismatch(es)\n'.format(
                varianttval))
        metamap['varianttval'] = varianttval
    else:
        metamap['varianttval'] = None

    # Compute Barcode Prefix t-value
    if not barcodeprefix is None:
        bpxtval = ci.infer_tvalue(
            elementlen=len(barcodeprefix),
            maximum=2)
        metamap['bpxtval'] = bpxtval
    else:
        metamap['bpxtval'] = None

    # Compute Barcode Suffix t-value
    if not barcodesuffix is None:
        bsxtval = ci.infer_tvalue(
            elementlen=len(barcodesuffix),
            maximum=2)
        metamap['bsxtval'] = bsxtval
    else:
        metamap['bsxtval'] = None

    # Compute Barcode Prefix t-value
    if not variantprefix is None:
        vpxtval = ci.infer_tvalue(
            elementlen=len(variantprefix),
            maximum=2)
        metamap['vpxtval'] = vpxtval
    else:
        metamap['vpxtval'] = None

    # Compute Barcode Suffix t-value
    if not variantsuffix is None:
        vsxtval = ci.infer_tvalue(
            elementlen=len(variantsuffix),
            maximum=2)
        metamap['vsxtval'] = vsxtval
    else:
        metamap['vsxtval'] = None

    # Compute Read Anchor Motif
    if not barcodeprefix is None:
         metamap['anchor'] = barcodeprefix
    else:
         metamap['anchor'] = barcodesuffix

    # Save and Delete metamap
    liner.send(' Writing Meta Map ...')
    opath = indexdir+'meta.map'
    ut.savedict(
        dobj=metamap,
        filepath=opath)
    indexqueue.put(opath)
    del metamap
    liner.send(' Writing    Meta Map  : Done\n')

    # Train Barcode Model
    liner.send(' Training Barcode Model ...')
    X,Y = ci.split_map(
        maptosplit=barcodedict)
    t0  = tt.time()
    barcodemodel = sy.Scry().train(
        X=X,
        Y=Y,
        t=barcodetval)

    # Save Barcode Model
    opath = indexdir+'barcode.model'
    ut.savemodel(
        mobj=barcodemodel,
        filepath=opath)
    indexqueue.put(opath)
    del barcodedict
    del barcodemodel
    liner.send(' Writing Barcode Model: Done\n')

    # Save and Delete variantdict
    liner.send(' Writing Variant Map ...')
    opath = indexdir+'variant.map'
    ut.savedict(
        dobj=variantdict,
        filepath=opath)
    indexqueue.put(opath)
    del variantdict
    liner.send(' Writing Variant Map  : Done\n')

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Indexing Completed
    indexqueue.put(None)

def index(
    barcodedata,
    barcodecol,
    indexfile,
    barcodeprefix=None,
    barcodesuffix=None,
    variantdata=None,
    variantcol=None,
    variantprefix=None,
    variantsuffix=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Indexing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Index]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass barcodedata Parsing and Validation
    (bardf,
    barcodedata_valid) = vp.get_parsed_indata_info(
        indata=barcodedata,
        indata_field='  Barcode Data  ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Sort bardf
    if barcodedata_valid:
        bardf.sort_index(
            inplace=True)

    # Full barcodecol Validation
    (barcodes,
    barcodecol_valid) = vp.get_parsed_column_info(
        col=barcodecol,
        df=bardf,
        col_field='  Barcode Column',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full indexfile Validation
    indexfile_valid = vp.get_outfile_validity(
        outfile=indexfile,
        outfile_suffix='.oligopool.index',
        outfile_field='    Index File  ',
        liner=liner)

    # Adjust indexfile Suffix
    indexfile = ut.get_adjusted_path(
        path=indexfile,
        suffix='.oligopool.index')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full barcodeprefix Validation
    (barcodeprefix,
    barcodeprefix_valid) = vp.get_parsed_column_info(
        col=barcodeprefix,
        df=bardf,
        col_field='  Barcode Prefix',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full barcodesuffix Validation
    (barcodesuffix,
    barcodesuffix_valid) = vp.get_parsed_column_info(
        col=barcodesuffix,
        df=bardf,
        col_field='  Barcode Suffix',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # First Pass indata Parsing and Validation
    (vardf,
    variantdata_valid) = vp.get_parsed_variantdata_info(
        variantdata=variantdata,
        variantdata_field='  Variant Data  ',
        required_fields=('ID',),
        bardf=bardf,
        barcodedata_valid=barcodedata_valid,
        liner=liner)

    # Full variantcol Validation
    (variants,
    variantcol_valid) = vp.get_parsed_column_info(
        col=variantcol,
        df=vardf,
        col_field='  Variant Column',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full variantprefix Validation
    (variantprefix,
    variantprefix_valid) = vp.get_parsed_column_info(
        col=variantprefix,
        df=vardf,
        col_field='  Variant Prefix',
        col_desc='Input from Column',
        col_type=0,
        adjcol=variantcol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full variantsuffix Validation
    (variantsuffix,
    variantsuffix_valid) = vp.get_parsed_column_info(
        col=variantsuffix,
        df=vardf,
        col_field='  Variant Suffix',
        col_desc='Input from Column',
        col_type=0,
        adjcol=variantcol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # First Pass Validation
    if not all([
        barcodedata_valid,
        barcodecol_valid,
        indexfile_valid,
        barcodeprefix_valid,
        barcodesuffix_valid,
        variantdata_valid,
        variantcol_valid,
        variantprefix_valid,
        variantsuffix_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Setup Warning Dictionary
    warns = {}

    # Parse Barcode Data Feasibility
    liner.send('\n[Step 1: Parsing Barcode Data]\n')

    # Parse bardf
    (parsestatus,
    IDdict,
    barcodedict,
    duplicates,
    barcodecount,
    barcodesuniq,
    barcodelen,
    barcodelenuniq,
    minbarcodelen,
    maxbarcodelen) = ci.get_parsed_barcodes(
        bardf=bardf,
        barcodes=barcodes,
        liner=liner)

    # Barcode Data infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'stepname': 'parsing-barcode-data',
            'vars'    : {
                      'duplicates': duplicates,
                    'barcodecount': barcodecount,
                  'barcodesunique': barcodesuniq,
                   'minbarcodelen': minbarcodelen,
                   'maxbarcodelen': maxbarcodelen,
                'barcodelenunique': barcodelenuniq},
            'warns'   : warns}

        # Return results
        liner.close()
        return stats

    # Parse Barcode Constant Feasibility
    liner.send('\n[Step 2: Parsing Barcode Constants]\n')

    # Parse barcodeprefix and barcodesuffix
    (parsestatus,
    constantsextracted,
    constantsuniq,
    longconstants,
    barcodeprefix,
    barcodesuffix,
    prefixuniq,
    suffixuniq,
    prefixlen,
    suffixlen) = ci.get_parsed_constants(
        prefixconstant=barcodeprefix,
        suffixconstant=barcodesuffix,
        liner=liner)

    # Barcode Constants infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'stepname': 'parsing-barcode-constants',
            'vars'    : {
                'constantsextracted': constantsextracted,
                   'constantsunique': constantsuniq,
                     'longconstants': longconstants,
                      'prefixunique': prefixuniq,
                      'suffixunique': suffixuniq,
                         'prefixlen': prefixlen,
                         'suffixlen': suffixlen},
            'warns'   : warns}

        # Return results
        liner.close()
        return stats

    # Do we have associated variants?
    if not variantdata is None:

        # Extract Variant Data
        liner.send('\n[Step 3: Extracting Variant Data]\n')

        # Update Step 3 Warning
        warns[3] = {
            'warncount': 0,
            'stepname' : 'parsing-variant-data',
            'vars': None}

        # Extract variantdata
        variantdict = ci.get_extracted_variants(
            variants=variants,
            expcount=barcodecount,
            IDdict=IDdict,
            warn=warns[3],
            liner=liner)

        # Parse Variant Constant Feasibility
        liner.send('\n[Step 4: Parsing Variant Constants]\n')

        # Parse variantprefix and variantsuffix
        (parsestatus,
        constantsextracted,
        constantsuniq,
        longconstants,
        variantprefix,
        variantsuffix,
        prefixuniq,
        suffixuniq,
        prefixlen,
        suffixlen) = ci.get_parsed_constants(
            prefixconstant=variantprefix,
            suffixconstant=variantsuffix,
            liner=liner)

        # Variant Constants infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
                'stepname': 'parsing-variant-constants',
                'vars'    : {
                    'constantsextracted': constantsextracted,
                    'constantsunique': constantsuniq,
                        'longconstants': longconstants,
                        'prefixunique': prefixuniq,
                        'suffixunique': suffixuniq,
                            'prefixlen': prefixlen,
                            'suffixlen': suffixlen},
                'warns'   : warns}

            # Return results
            liner.close()
            return stats

    else:
        variantdict   = {}
        variantprefix = None
        variantsuffix = None

    # Setup Workspace
    (indexfile,
    indexdir) = ut.setup_workspace(
        outfile=indexfile,
        outfile_suffix='.oligopool.index')

    # Schedule indexfile deletion
    ifdeletion = ae.register(
        ut.remove_file,
        indexfile)

    # Prepared Objects Queue
    indexqueue = mp.SimpleQueue()

    # Launching Indexing
    liner.send('\n[Step 5: Computing Index]\n')

    # Define Index Stats
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 5,
        'stepname': 'computing-index',
        'vars'    : None,
        'warns'   : warns}

    # Compute Index Objects
    index_engine(
        IDdict=IDdict,
        barcodedict=barcodedict,
        barcodename=barcodecol,
        barcodecount=barcodecount,
        barcodelen=barcodelen,
        barcodeprefix=barcodeprefix,
        barcodesuffix=barcodesuffix,
        variantdict=variantdict,
        variantprefix=variantprefix,
        variantsuffix=variantsuffix,
        indexdir=indexdir,
        indexqueue=indexqueue,
        liner=liner)

    # Archive All Indexed Objects
    ut.archive(
        objqueue=indexqueue,
        arcfile=indexfile,
        prodcount=1,
        prodactive=0,
        liner=liner)

    # Indexing Status?
    indexstatus = 'Successful'

    # Indexing Stats
    liner.send('\n[Indexing Stats]\n')

    liner.send(
        ' Indexing Status : {}\n'.format(
            indexstatus))
    liner.send(
        '     Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Remove Workspace
    ut.remove_directory(
        dirpath=indexdir)

    # Unschedule packfile deletion
    if indexstatus == 'Successful':
        ae.unregister(ifdeletion)

    # Close Liner
    liner.close()

    return stats