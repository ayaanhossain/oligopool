import time  as tt

import collections as cx
import numpy       as np
import atexit      as ae

import utils       as ut
import valparse    as vp
import corebarcode as cb


def barcode_engine(
    barcodelength,
    minhdist,
    exmotifs,
    leftcontext,
    rightcontext,
    edgeeffectlength,
    targetcount,
    stats,
    liner):
    '''
    Return barcodes fulfilling constraints.
    Internal use only.

    :: barcodelength
       type - integer
       desc - required barcode length
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: leftcontext
       type - list / None
       desc - list of sequence to the
              left of barcodes,
              None otherwise
    :: rightcontext
       type - list / None
       desc - list of sequences to the
              right of barcodes,
              None otherwise
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    :: targetcount
       type - integer
       desc - required number of barcodes
              to be designed
    :: stats
       type - dict
       desc - barcode design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()              # Start Timer
    store = np.zeros(           # Store Encoded Barcodes
        (targetcount, barcodelength),
        dtype=np.float64)
    codes = []                  # Store Decoded Barcodes
    assignmentarray = []        # Assignment Array

    # Context Setup
    contextarray  = cx.deque(range(targetcount))   # Context Array
    (leftcontexttype,
    leftselector)  = cb.get_context_type_selector( # Left  Context Selector
        context=leftcontext)
    (rightcontexttype,
    rightselector) = cb.get_context_type_selector( # Right Context Selector
        context=rightcontext)

    # Setup Jumper and Generator
    (jtp,
    jumper)  = cb.get_jumper(
        barcodelength=barcodelength)
    barcodes = cb.stream_barcodes(
        barcodelength=barcodelength,
        jumper=jumper)
    barcode    = None # Current Barcode
    acccodeseq = None # Last Successful Barcode Sequence

    # Infinite Jumper Failure
    prob  = ut.get_prob(    # Probability of Success
        success=1,
        trials=4 ** (
            barcodelength - minhdist))
    trial = ut.get_trials(  # Trials Required
        prob=prob)
    sscnt = 1     # Success Count
    flcnt = 0     # Failure Count
    excnt = trial # Trial   Count

    # Verbage Setup
    verbage_reach  = 0
    verbage_target = ut.get_sample(
        value=targetcount,
        lf=0.080,
        uf=0.120)
    plen = ut.get_printlen(
        value=targetcount)

    # Build Barcodes
    while True:

        # Sample a Barcode in Space
        barcode = next(barcodes)

        # Space Exhausted?
        if barcode is None:

            # Final Update
            cb.show_update(
                count=stats['vars']['barcodecount'],
                plen=plen,
                barcodeseq=acccodeseq,
                optstatus=True,
                optstate=0,
                inittime=t0,
                terminal=True,
                liner=liner)

            # No Solution
            return (None,
                None,
                stats)

        # Update Barcode Store
        store[stats['vars']['barcodecount'], :] = barcode

        # Decode Barcode Sequence
        barcodeseq = cb.get_barcodeseq(
            barcode=barcode)

        # Check Feasibility
        (optstatus,
        optstate,
        aidx,
        failcount) = cb.is_barcode_feasible(
            store=store,
            count=stats['vars']['barcodecount'],
            barcodeseq=barcodeseq,
            barcodelength=barcodelength,
            minhdist=minhdist,
            exmotifs=exmotifs,
            leftcontexttype=leftcontexttype,
            leftselector=leftselector,
            rightcontexttype=rightcontexttype,
            rightselector=rightselector,
            edgeeffectlength=edgeeffectlength,
            contextarray=contextarray)

        # Inifinite Jumper Book-keeping Update
        if jtp == 2:
            excnt += 1

        # Accept Sample into Store?
        if optstatus:

            # Update Store Fill Count
            stats['vars']['barcodecount'] += 1

            # Show Update on Success
            cb.show_update(
                count=stats['vars']['barcodecount'],
                plen=plen,
                barcodeseq=barcodeseq,
                optstatus=optstatus,
                optstate=optstate,
                inittime=t0,
                terminal=False,
                liner=liner)

            # Record Assignment Index
            assignmentarray.append(aidx)

            # Update Codes
            codes.append(barcodeseq)

            # Update Last Accounted Barcode
            acccodeseq = barcodeseq

            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                sscnt += 1
                prob  = ut.get_prob(
                    success=sscnt,
                    trials=excnt)
                trial = ut.get_trials(
                    prob=prob)
                flcnt = 0

        # Rejected from Store
        else:

            # Update Fail Counts
            if optstate == 1:
                stats['vars']['distancefail'] += 1
            if optstate == 2:
                stats['vars']['motiffail']    += 1
            if optstate == 3:
                stats['vars']['edgefail']     += 1

            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                flcnt += 1

            # Record Context Failure Motifs
            if not failcount is None:
                stats['vars']['motifcounter'] += failcount

        # Verbage Book-keeping
        if verbage_reach >= verbage_target:
            cb.show_update(
                count=stats['vars']['barcodecount'],
                plen=plen,
                barcodeseq=barcodeseq,
                optstatus=optstatus,
                optstate=optstate,
                inittime=t0,
                terminal=False,
                liner=liner)
            verbage_reach = -1
        verbage_reach += 1

        # Target Reached?
        if stats['vars']['barcodecount'] == targetcount:

            # Construct the Sorted Barcodes
            codes = [code for aidx,code in sorted(
                zip(assignmentarray, codes))]

            # Update Design Status
            stats['status'] = True
            stats['basis']  = 'solved'

            # Final Update
            cb.show_update(
                count=stats['vars']['barcodecount'],
                plen=plen,
                barcodeseq=codes[-1],
                optstatus=True,
                optstate=0,
                inittime=t0,
                terminal=True,
                liner=liner)

            # Return Solution
            return (codes,
                store,
                stats)

        # Trials Exhausted for Inifinite Jumper?
        if jtp == 2:
            if flcnt == trial:

                # Final Update
                cb.show_update(
                    count=stats['vars']['barcodecount'],
                    plen=plen,
                    barcodeseq=acccodeseq,
                    optstatus=True,
                    optstate=0,
                    inittime=t0,
                    terminal=True,
                    liner=liner)

                # No Solution
                return (None,
                    None,
                    stats)

def barcode(
    indata,
    barcodelength,
    minhdist,
    barcodecol,
    outfile=None,
    exmotifs=None,
    leftcontext=None,
    rightcontext=None,
    verbose=True):
    '''
    The barcode function designs constrained barcodes, such that
    the minimum hamming distance between every pair of barcodes
    in the set is guaranteed to be greater than or equal to the
    specified minimum hamming distance. Additionally, the set of
    barcodes are free of all specified motifs, even when placed
    between the left and right context sequences. The generated
    DataFrame containing designed barcodes is returned to user
    and optionally written out to <outfile> (CSV) if specified.

    :: indata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas DataFrame storing
              annotated oligopool variants and their parts
    :: barcodelength
       type - integer
       desc - the length of the barcodes to be designed,
              must be 4 or greater
    :: minhdist
       type - integer
       desc - minimum required hamming distance between every
              pair of barcodes in the designed set, must be 1
              or greater
    :: barcodecol
       type - string
       desc - the name of the column to store designed barcodes
    :: outfile
       type - string
       desc - filename to save updated DataFrame with barcodes
              (suffix='.oligopool.barcode.csv')
    :: exmotifs
       type - iterable / string / pd.DataFrame / None
       desc - iterable of DNA string motifs to be excluded
              within and at the edges of the barcodes when
              placed around context sequences; optionally,
              this can be a path to a CSV file containing
              uniquely identified excluded motifs or an
              equivalent pandas DataFrame
              (default=None)
    :: leftcontext
       type - string / None
       desc - column name containing the DNA sequences placed
              left of the barcode sequence in 5' to 3' direction
              of the forward strand
              (default=None)
    :: rightcontext
       type - string / None
       desc - column name containing the DNA sequences placed
              right of the barcode sequence in 5' to 3' direction
              of the forward strand
              (default=None)
    :: verbose
       type - boolean
       desc - if True will log updates to stdout
              (default=True)

    Output: A file named <outfile> with '.oligoool.barcode.csv'
            suffix if specified; otherwise a pandas DataFrame is
            returned.

    Note 1. Specified <indata> must contain a column named 'ID',
            that uniquely identifies variants in the pool.

    Note 2. If <exmotifs> points to a  CSV file or DataFrame,
            it must contain both an 'ID' column and an 'Exmotifs'
            column, with 'Exmotifs' containing excluded motifs.

    Note 3. All column names in <indata> must be unique, without
            <barcodename> as a pre-existing column name.

    Note 4. All rows and columns in <indata> must be non-empty,
            i.e. none of the cells must be empty.

    Note 5. The columns <leftcontext> and <rightcontext> must be
            adjacent to each other and in order.

    Note 6. All values in <indata> except those in 'ID' must be
            DNA strings.
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Barcode]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_data_info(
        data=indata,
        data_field='    Input Data    ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full barcodelength Validation
    barcodelength_valid = vp.get_numeric_validity(
        numeric=barcodelength,
        numeric_field='  Barcode Length  ',
        numeric_pre_desc=' ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full minhdist Validation
    minhdist_valid = vp.get_numeric_validity(
        numeric=minhdist,
        numeric_field='  Hamming Distance',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' Mismatch(es) per Barcode Pair',
        minval=1,
        maxval=barcodelength if barcodelength_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full outcol Validation
    barcodecol_valid = vp.get_parsed_column_info(
        col=barcodecol,
        df=indf,
        col_field='  Barcode Column  ',
        col_desc='Output in Column',
        col_type=1,
        adjcol=None,
        adjval=None,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.barcode.csv',
        outdf_field='   Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.barcode.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field=' Excluded Motifs  ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotifs,',
        required=False,
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='     Left Context ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=rightcontextname,
        adjval=+1,
        liner=liner)

    # Full leftcontext Parsing and Validation
    (rightcontext,
    rightcontext_valid) = vp.get_parsed_column_info(
        col=rightcontext,
        df=indf,
        col_field='    Right Context ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        barcodelength_valid,
        minhdist_valid,
        barcodecol_valid,
        outfile_valid,
        exmotifs_valid,
        leftcontext_valid,
        rightcontext_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    barcodelength = round(barcodelength)
    minhdist      = round(minhdist)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Barcode Design Book-keeping
    outdf = None
    stats = None

    # Parse Barcode Length Feasibility
    liner.send('\n[Parsing Barcode Length]\n')

    # Parse barcodelength
    (parsestatus,
    designspace,
    targetcount) = cb.get_parsed_barcode_length(
        barcodelength=barcodelength,
        indf=indf,
        liner=liner)

    # barcodelength infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 1,
            'vars'  : {
                'barcodelength': barcodelength,
                        'designspace': designspace,
                        'targetcount': targetcount}}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Parsing Excluded Motifs]\n')

        # Parse exmotifs
        (parsestatus,
        exmotifs,
        problens) = ut.get_parsed_exmotifs(
            exmotifs=exmotifs,
            typer=tuple,
            element='Barcode',
            liner=liner)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status': False,
                'basis' : 'infeasible',
                'step'  : 2,
                'vars'  : {
                     'problens': problens,
                    'probcount': tuple(list(
                        4**pl for pl in problens))}}

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            exmotifs=exmotifs)

    # Extract Left and Right Context
    if not exmotifs is None and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Show update
        liner.send('\n[Extracting Context Sequences]\n')

        # Extract Both Contexts
        (leftcontext,
        rightcontext) = ut.get_extracted_context(
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            edgeeffectlength=edgeeffectlength,
            reduce=False,
            liner=liner)
    else:
        (leftcontext,
        rightcontext) = None, None

    # Launching Barcode Design
    liner.send('\n[Computing Barcodes]\n')

    # Define Barcode Design Stats
    stats = {
        'status': False,
         'basis': 'unsolved',
          'step': 4,
          'vars': {
               'targetcount': targetcount,   # Required Number of Barcodes
              'barcodecount': 0,             # Barcode Design Count
              'distancefail': 0,             # Hamming Distance Fail Count
                 'motiffail': 0,             # Motif Elimination Fail Count
                  'edgefail': 0,             # Edge Effect Fail Count
            'distancedistro': None,          # Hamming Distance Distribution
              'motifcounter': cx.Counter()}} # Motif Encounter Counter

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Barcodes
    (codes,
    store,
    stats) = barcode_engine(
        barcodelength=barcodelength,
        minhdist=minhdist,
        exmotifs=exmotifs,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        edgeeffectlength=edgeeffectlength,
        targetcount=targetcount,
        stats=stats,
        liner=liner)

    # Compute Hamming Distribution
    if not store is None:
        liner.send('\n[Computing Distance Distribution]\n')
        stats['vars']['distancedistro'] = cb.get_distro(
            store=store,
            codes=codes,
            liner=liner)

    # Barcode Status
    if stats['status']:
        barcodestatus = 'Successful'
    else:
        barcodestatus = 'Failed'

    # Insert codes into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=codes,
            outcol=barcodecol)

        # Prepare outdf
        outdf = indf

        # Write indf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Barcoding Statistics
    liner.send('\n[Barcode Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'targetcount',
            'barcodecount')))

    liner.send(
        '   Design Status   : {}\n'.format(
            barcodestatus))
    liner.send(
        '   Target Count    : {:{},d} Barcode(s)\n'.format(
            stats['vars']['targetcount'],
            plen))
    liner.send(
        ' Barcodes Count    : {:{},d} Barcode(s) ({:6.2f} %)\n'.format(
            stats['vars']['barcodecount'],
            plen,
            ut.safediv(
                A=stats['vars']['barcodecount'] * 100.,
                B=targetcount)))

    # Success Relevant Stats
    if stats['status']:
        if stats['vars']['distancedistro']:

            dlen = ut.get_printlen(
                value=max(stats['vars']['distancedistro'].keys()))

            clen = ut.get_printlen(
                value=max(stats['vars']['distancedistro'].values()))

            liner.send('   Pair-wise Distance Distribution\n')

            for distance,count in stats['vars']['distancedistro'].most_common():
                liner.send(
                    '     - {:{},d} Barcode(s) w/ Distance ≥ {:{},d} Mismatches\n'.format(
                        count,
                        clen,
                        distance,
                        dlen))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'distancefail',
            'motiffail',
            'edgefail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['distancefail'] + \
                            stats['vars']['motiffail']    + \
                            stats['vars']['edgefail']
        liner.send(
            ' Distance Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['distancefail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['distancefail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '    Motif Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['motiffail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['motiffail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '     Edge Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['edgefail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['edgefail'] * 100.,
                    B=total_conflicts)))

        # Enumerate Motif-wise Fail Counts
        if stats['vars']['motifcounter']:

            qlen = max(len(motif) \
                for motif in stats['vars']['motifcounter'].keys()) + 2

            sntn, vlen = ut.get_notelen(
                printlen=ut.get_printlen(
                    value=max(
                        stats['vars']['motifcounter'].values())))

            liner.send('   Motif-wise Conflict Distribution\n')

            for motif,count in stats['vars']['motifcounter'].most_common():
                motif = '\'{}\''.format(motif)
                liner.send(
                    '     - Motif {:>{}} Triggered {:{},{}} Event(s)\n'.format(
                        motif, qlen, count, vlen, sntn))

    # Show Time Elapsed
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Unschedule outfile deletion
    if barcodestatus:
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
