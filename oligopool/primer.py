import time  as tt

import collections as cx
import atexit      as ae

import nrpcalc    as nr

import utils      as ut
import valparse   as vp
import coreprimer as cp


def primer_objectives(
    primer,
    primerlen,
    primertype,
    mintmelt,
    maxtmelt,
    maxreplen,
    oligorepeats,
    pairedprimer,
    pairedrepeats,
    exmotifs,
    prefixforbidden,
    suffixforbidden,
    background,
    stats,
    liner):
    '''
    TBW
    '''

    # Objective 1: Background Non-Repetitiveness
    obj1, traceloc = cp.is_background_feasible(
        primer=primer,
        background=background)

    # Objective 1 Failed
    if not obj1:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Background Repeat'.format(
                primer))

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 2: Oligopool Non-Repetitiveness
    obj2, traceloc = cp.is_oligopool_feasible(
        primer=primer,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats)

    # Objective 2 Failed
    if not obj2:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Oligopool Repeat'.format(
                primer))

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 3: Paired Primer Non-Repetitiveness
    obj3, traceloc = cp.is_paired_feasible(
        primer=primer,
        pairedrepeats=pairedrepeats)

    # Objective 3 Failed
    if not obj3:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Paired Repeat'.format(
                primer))

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 4: Melting Temeperature Bounded
    obj4, failtype = cp.is_tmelt_feasible(
        primer=primer,
        primerlen=primerlen,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt)

    # Objective 4 Failed
    if not obj4:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Tm Infeasibility'.format(
                primer))

        # Compute Traceback
        traceloc = cp.get_tmelt_traceback(
            primer=primer,
            failtype=failtype)

        # Update Stats
        stats['vars']['Tmfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 5: Motif Embedding
    obj5, motif = cp.is_motif_feasible(
        primer=primer,
        exmotifs=exmotifs)

    # Objective 5 Failed
    if not obj5:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Motif Infeasibility'.format(
                primer))

        # Update Stats
        stats['vars']['motiffail'] += 1
        stats['vars']['motifcounter'][motif] += 1

        # Return Traceback
        return False, max(
            0,
            len(primer)-1)

    # Objective 6: Edge Feasibility (Edge-Effects)
    obj6, traceloc = cp.is_edge_feasible(
        primer=primer,
        primerlen=primerlen,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden)

    # Objective 6 Failed
    if not obj6:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Edge Infeasibility'.format(
                primer))

        # Update Stats
        stats['vars']['edgefail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 7: Homodimer Feasibility
    obj7, traceloc = cp.is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        primerspan=None,
        pairedprimer=None,
        pairedspan=None,
        dimertype=0)

    # Objective 7 Failed
    if not obj7:

        # Show Update
        liner.send(
            ' Candidate: Primer {} Rejected due to Homodimer Infeasibility'.format(
                primer))

        # Update Stats
        stats['vars']['homodimerfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 8: Heterodimer Feasibility
    obj8, traceloc = cp.is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        primerspan=None,
        pairedprimer=pairedprimer,
        pairedspan=None,
        dimertype=1)

    # Objective 8 Failed
    if not obj8:

        # Show Update
        liner.send(
            ' Candidate: {} Rejected due to Heterodimer Infeasibility'.format(
                primer))

        # Update Stats
        stats['vars']['heterodimerfail'] += 1

        # Return Traceback
        return False, traceloc

    # Show Update
    if len(primer) < primerlen:
        liner.send(
            ' Candidate: Primer {} is Partially Accepted'.format(
                primer))
    else:
        liner.send(
            ' Candidate: Primer {} is Accepted'.format(
                primer))

    # All Objectives OK!
    return True

def primer_engine(
    primerseq,
    primertype,
    mintmelt,
    maxtmelt,
    maxreplen,
    oligorepeats,
    pairedprimer,
    pairedrepeats,
    exmotifs,
    prefixforbidden,
    suffixforbidden,
    background,
    stats,
    liner):
    '''
    TBW
    '''

    # Book-keeping
    primerstruct = 'x'*len(primerseq)    # Primer Structure Constraint
    t0 = tt.time()

    # Update Orientation
    # Since Forward Primer Design,
    # interprest Paired Primer as
    # Reverse Primer specified in
    # terns of Forward Strand
    if primertype == 0:
        if not pairedprimer is None:
            pairedprimer = ut.get_revcomp(
                seq=pairedprimer)

    # Optimize exmotifs
    exmotifs = ut.get_grouped_sequences(
        sequences=exmotifs)

    # Define Objective Function
    objectivefunction = lambda primer: primer_objectives(
        primer=primer,
        primerlen=len(primerseq),
        primertype=primertype,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=pairedprimer,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden,
        background=background,
        stats=stats,
        liner=liner)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Design Primer via Maker
    primer = maker.nrp_maker(
        homology=6,
        seq_constr=primerseq,
        struct_constr=primerstruct,
        target_size=1,
        background=None,
        struct_type='both',
        synth_opt=True,
        local_model_fn=objectivefunction,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Show Time Elapsed
    liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
        tt.time() - t0))

    # Check Status and Return Solution
    if len(primer) > 0:

        # Extract Primer
        primer = primer[0]

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

        # Update Tm
        stats['vars']['primerTm'] = ut.get_tmelt(
            seq=primer)

        # Update GC Percentage
        stats['vars']['primerGC'] = (primer.count('G') + \
                                     primer.count('C')) \
                                        / (len(primer) * 0.01)

        # Update Hairpin Free Energy
        stats['vars']['hairpinMFE'] = cp.folder.evaluate_mfe(
            seq=primer,
            dg=True)[-1]

        # Update Heterodimer Free Energy
        stats['vars']['homodimerMFE'] = cp.folder.evaluate_mfe_dimer(
            seq1=primer,
            seq2=primer)[-1]

        # Update Homodimer Free Energy
        if not pairedprimer is None:
            stats['vars']['heterodimerMFE'] = cp.folder.evaluate_mfe_dimer(
                seq1=primer,
                seq2=pairedprimer)[-1]

        # Return Results
        return (primer, stats)

    # Design Unsuccessful
    else:
        return (None, stats)

def primer(
    indata,
    primerseq,
    primertype,
    mintmelt,
    maxtmelt,
    maxreplen,
    primercol,
    outfile=None,
    pairedcol=None,
    exmotifs=None,
    leftcontext=None,
    rightcontext=None,
    background=None,
    verbose=True):
    '''
    The primer function designs constrained primers, at desired
    melting temperature, with desired non-repetitiveness, that
    works optimally for all variants in the oligopool, without
    any excluded motifs inside or introducing new ones at the
    primer edges. Additional constraints are enforced to ensure
    compatibility with a paired primer. The generated DataFrame
    containing designed primer is returned, and optionally also
    written out to <outfile> (CSV) if specified.

    :: indata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas DataFrame storing
              annotated oligopool variants and their parts
    :: primerseq
       type - integer
       desc - an IUPAC degenerate primer sequence constraint
    :: primertype
       type - integer
       desc - primer design type identifier
              0 = a forward primer is designed
              1 = a reverse primer is designed
    :: mintmelt
       type - float
       desc - minimum allowed primer melting temperature in
              degree celsium, must be 25 or greater
    :: maxtmelt
       type - float
       desc - maximum allowed primer melting temperature in
              degree celsium, must be 95 or lesser
    :: maxreplen
       type - integer
       desc - maximum shared repeat length between the primers
              and flanking regions, must be between 6 and 20
    :: primercol
       type - string
       desc - the name of the column to store designed primer
    :: outfile
       type - string
       desc - filename to save updated DataFrame with primers
              (suffix='.oligopool.primer.csv')
    ::
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Primer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_data_info(
        data=indata,
        data_field='      Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # First Pass primerseq Validation
    primerseq_valid = vp.get_primerseq_validity(
        primerseq=primerseq,
        primerseq_field='     Primer Sequence   ',
        liner=liner)

    # Full primertype Validation
    primertype_valid = vp.get_categorical_validity(
        category=primertype,
        category_field='     Primer Type       ',
        category_pre_desc=' ',
        category_post_desc=' Primer Design',
        category_dict={
            0: 'Forward',
            1: 'Reverse'},
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    (mintmelt,
    maxtmelt,
    tmelt_valid) = vp.get_parsed_range_info(
        minval=mintmelt,
        maxval=maxtmelt,
        range_field='    Melting Temperature',
        range_unit='°C',
        range_min=25,
        range_max=95,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length     ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=6,
        maxval=20,
        precheck=False,
        liner=liner)

    # Full primercol Validation
    primercol_valid = vp.get_parsed_column_info(
        col=primercol,
        df=indf,
        col_field='     Primer Column     ',
        col_desc='Output in Column',
        col_type=1,
        adjcol=None,
        adjval=None,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.primer.csv',
        outdf_field='     Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.primer.csv')

    # # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full pairedprimer Validation
    (pairedprimer,
    pairedprimer_valid) = vp.get_constantcol_validity(
        constantcol=pairedcol,
        constantcol_field='     Paired Primer     ',
        df=indf,
        liner=liner)

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field='   Excluded Motifs     ',
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
        col_field='       Left Context    ',
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
        col_field='      Right Context    ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        liner=liner)

    # Full background Parsing and Validation
    (background,
    background_valid) = vp.get_parsed_background(
        background=background,
        background_field=' Background Database   ',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        primerseq_valid,
        primertype_valid,
        tmelt_valid,
        maxreplen_valid,
        primercol_valid,
        outfile_valid,
        pairedprimer_valid,
        exmotifs_valid,
        leftcontext_valid,
        rightcontext_valid,
        background_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust maxreplen
    maxreplen = round(maxreplen)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Barcode Design Book-keeping
    outdf = None
    stats = None

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
            element='Primer',
            liner=liner)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status': False,
                'basis' : 'infeasible',
                'step'  : 1,
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

    # Parsing Sequence Constraint Feasibility
    liner.send('\n[Parsing Primer Sequence]\n')

    # Full primerseq Validation
    (parsestatus,
    designspace,
    excludedmotifs,
    internalrepeats,
    palindromes,
    pairedrepeats) = cp.get_parsed_sequence_constraint(
        primerseq=primerseq,
        primertype=primertype,
        exmotifs=exmotifs,
        pairedprimer=pairedprimer,
        liner=liner)

    # primerseq infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 2,
            'vars'  : {
                    'designspace': designspace,
                 'excludedmotifs': excludedmotifs,
                'internalrepeats': internalrepeats,
                    'palindromes': palindromes,
                  'pairedrepeats': pairedrepeats}}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define pairedrepeats
    if not pairedprimer is None:
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=pairedprimer,
            k=6))
    else:
        pairedrepeats = None

    # Parse Melting Temperature
    liner.send('\n[Parsing Melting Temperature]\n')

    # Full mintmelt and maxtmelt Validation
    (parsestatus,
    estimatedminTm,
    estimatedmaxTm,
    higherminTm,
    lowermaxTm,
    mintmelt,
    maxtmelt) = cp.get_parsed_tmelt_constraint(
        primerseq=primerseq,
        pairedprimer=pairedprimer,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        liner=liner)

    # primerseq infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 3,
            'vars'  : {
                'estimatedminTm': estimatedminTm,
                'estimatedmaxTm': estimatedmaxTm,
                   'higherminTm': higherminTm,
                    'lowermaxTm': lowermaxTm}}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Edge Effects
    if not exmotifs is None and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Show update
        liner.send('\n[Extracting Context Sequences]\n')

        # Extract Primer Contexts
        (leftcontext,
        rightcontext) = ut.get_extracted_context(
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            edgeeffectlength=edgeeffectlength,
            reduce=False,
            liner=liner)

        # Show update
        liner.send('\n[Parsing Edge Effects]\n')

        # Parse primerseq, exmotifs,
        # leftcontext and rightcontext
        (parsestatus,
        probseqprefixes,
        probprefixlens,
        probseqsuffixes,
        probsuffixlens,
        prefixforbidden,
        suffixforbidden) = cp.get_parsed_edgeeffects(
            primerseq=primerseq,
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            exmotifs=exmotifs,
            liner=liner)

        # Edge Effect infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status': False,
                'basis' : 'infeasible',
                'step'  : 5,
                'vars'  : {
                    'probseqprefixes': probseqprefixes,
                     'probprefixlens': probprefixlens,
                    'probseqsuffixes': probseqsuffixes,
                     'probsuffixlens': probsuffixlens}}

            # Return results
            liner.close()
            return (outdf, stats)

    else:
        (prefixforbidden,
        suffixforbidden) = None, None

    # Extract Oligopool Repeats
    liner.send('\n[Extractng Oligopool Repeats]\n')

    # Extract Repeats from indf
    (parsestatus,
    kmerspace,
    fillcount,
    leftcount,
    oligorepeats) = ut.get_parsed_oligo_repeats(
        df=indf,
        maxreplen=maxreplen,
        element='Primer',
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 6,
            'vars'  : {
                'kmerspace': kmerspace,
                'fillcount': fillcount,
                'leftcount': leftcount}}

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Barcode Design
    liner.send('\n[Computing Primer]\n')

    # Define Barcode Design Stats
    stats = {
        'status': False,
         'basis': 'unsolved',
          'step': 7,
          'vars': {
                   'primerTm': None,          # Primer Melting Temperature
                   'primerGC': None,          # Primer GC Content
                 'hairpinMFE': None,          # Primer Hairpin Free Energy
               'homodimerMFE': None,          # Homodimer Free Energy
             'heterodimerMFE': None,          # Heterodimer Free Energy
                     'Tmfail': 0,             # Melting Temperature Fail Count
                 'repeatfail': 0,             # Repeat Fail Count
              'homodimerfail': 0,             # Homodimer Fail Count
            'heterodimerfail': 0,             # Heterodimer Fail Count
                  'motiffail': 0,             # Motif Elimination Fail Count
                   'edgefail': 0,             # Edge Effect Fail Count
               'motifcounter': cx.Counter()}} # Motif Encounter Counter

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Primer
    (primer,
    stats) = primer_engine(
        primerseq=primerseq,
        primertype=primertype,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=pairedprimer,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden,
        background=background,
        stats=stats,
        liner=liner)

    # Primer Status
    if stats['status']:
        primerstatus = 'Successful'
    else:
        primerstatus = 'Failed'

    # Insert primer into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=primer,
            outcol=primercol)

        # Prepare outdf
        outdf = indf

        # Write indf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Barcoding Statistics
    liner.send('\n[Primer Design Statistics]\n')

    liner.send(
        '      Design Status     : {}\n'.format(
            primerstatus))

    # Success Relevant Stats
    if stats['status']:

        liner.send(
            '     Melting Temperature: {:6.2f} °C\n'.format(
                stats['vars']['primerTm']))
        liner.send(
            '          GC Content    : {:6.2f} %\n'.format(
                stats['vars']['primerGC']))
        liner.send(
            '     Hairpin MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['hairpinMFE']))
        liner.send(
            '   Homodimer MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['homodimerMFE']))

        if not pairedprimer is None:
            liner.send(
                ' Heterodimer MFE        : {:6.2f} kcal/mol\n'.format(
                    stats['vars']['heterodimerMFE']))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'Tmfail',
            'repeatfail',
            'homodimerfail',
            'heterodimerfail',
            'motiffail',
            'edgefail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['Tmfail']          + \
                          stats['vars']['repeatfail']      + \
                          stats['vars']['homodimerfail']   + \
                          stats['vars']['heterodimerfail'] + \
                          stats['vars']['motiffail']       + \
                          stats['vars']['edgefail']

        liner.send(
            ' Melt. Temp. Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['Tmfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['Tmfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '      Repeat Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeatfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeatfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '   Homodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['homodimerfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['homodimerfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            ' Heterodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['heterodimerfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['heterodimerfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '       Motif Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['motiffail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['motiffail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '        Edge Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
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
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))

    # Unschedule outfile deletion
    if primerstatus:
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
