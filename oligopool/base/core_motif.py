import time  as tt

import collections as cx

import numpy   as np
import nrpcalc as nr

from . import utils as ut


# Parser and Setup Functions

def get_parsed_sequence_constraint(
    motifseq,
    exmotifs,
    warn,
    liner):
    '''
    Check motif sequence feasibility.
    Internal use only.

    :: motifseq
       type - string
       desc - motif sequence constraint
    :: exmotifs
       type - deque / None
       desc - deque of all motifs
              to be excluded
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()
    optrequired  = True # Optimization Required
    exmotifindex = None # No Conflicts
    homology     = 6    # Initially, for Maker

    # Design Space Analysis
    liner.send(' Computing Design Space ...')

    dspace = 1
    for nt in motifseq:
        dspace *= len(ut.ddna_space[nt])
    sntn, plen = ut.get_notelen(
        printlen=ut.get_printlen(
            value=dspace))

    # Show Update
    if dspace == 1:
        optrequired = False # No Optimization .. Motif Constant
        liner.send(
            ' Design Space: 1 Possible Motif(s)\n')
    else:
        liner.send(
            ' Design Space: {:{},{}} Possible Motif(s)\n'.format(
                dspace,
                plen,
                sntn))

    # Exmotifs Analysis
    if dspace > 1 and not exmotifs is None:
        liner.send(' Computing Motif Conflicts ...')

        motif_ok, excludedmotifs = ut.get_exmotif_conflict(
            seq=motifseq,
            seqlen=len(motifseq),
            exmotifs=exmotifs,
            partial=False,
            checkall=True)

        # Show Update
        if not motif_ok:
            optrequired = True # Partial Motif Conflict Optimization

            # Update Warning Entry
            warn['vars'] = {'exmotifembedded': set()}

            # Compute Embedded Motif Indices
            # to be Ignored Downstream
            exmotlocdict = ut.get_exmotif_conflict_index(
                seq=motifseq,
                conflicts=excludedmotifs)
            exmotifindex = set()
            for exmotif in exmotlocdict:
                for loc in exmotlocdict[exmotif]:
                    exmotifindex.add(loc+len(exmotif))

            # Show Updates
            liner.send(
                ' Found {:,} Excluded Motif(s)\n'.format(
                    len(excludedmotifs)))

            # Record Warnings
            warn['warncount'] = len(excludedmotifs)
            warn['vars']['exmotifembedded'].update(excludedmotifs)
            homology = max(homology,
                           max(map(len, excludedmotifs)) + 1)

            # Show Excluded Motifs
            plen = max(map(len, excludedmotifs)) + 2
            for motif in excludedmotifs:
                motif = '\'{}\''.format(motif)
                liner.send(
                    '   - Excluded Motif {:>{}} Present [WARNING] (Excluded Motif Embedded)\n'.format(
                        motif,
                        plen))
        else:
            liner.send(
                ' Found 0 Excluded Motif(s)\n')

    # No Exmotifs Exist
    else:
        optrequired = False # No Optimization .. No Exmotifs

    # Region Analysis
    liner.send(' Computing Constant Regions ...')

    # Compute Constant Regions
    regions = ut.get_constant_regions(
        seqconstr=motifseq)

    # Finalize homology
    if regions:
        homology = max(homology,
                       max(map(len, regions)) + 1)

    # Show Time Elapsed
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Show Verdict
    if optrequired == 0:
        liner.send(
            ' Verdict: Motif Design is Constant\n')
    else:
        if exmotifindex:
            liner.send(
                ' Verdict: Motif Design with Embedded Excluded Motif(s)\n'
            )
        else:
            liner.send(
                ' Verdict: Motif Design Possibly Feasible\n')

    # Return Results
    return (optrequired, homology, exmotifindex)

def get_parsed_edgeeffects(
    motifseq,
    leftcontext,
    rightcontext,
    leftpartition,
    rightpartition,
    exmotifs,
    warn,
    liner):
    '''
    Cluster left and right context sequences
    and record forbidden prefix and suffixes.
    Internal use only.

    :: motifseq
       type - string
       desc - motif sequence constraint
    :: leftcontext
       type - tuple
       desc - tuple of all left context
              sequences
    :: rightcontext
       type - tuple
       desc - tuple of all right context
              sequences
    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    return ut.get_parsed_edgeeffects(
        sequence=motifseq,
        element='Motif',
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        leftpartition=leftpartition,
        rightpartition=rightpartition,
        exmotifs=exmotifs,
        merge=False,
        warn=warn,
        liner=liner)

def get_extracted_spacerlen(
    indf,
    oligolimit,
    liner):
    '''
    Extract spacer lengths based on existing
    variant length and maximum oligo length.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame containing
              all designed variants
    :: oligolimit
       type - integer
       desc - maximum allowed oligo length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Time-keeping
    t0 = tt.time()

    # Compute Variant Lengths
    liner.send(' Parsing Variant Lengths ...')

    variantlens = ut.get_variantlens(indf)

    plen = ut.get_printlen(
        value=max(len(str(index)) for index in indf.index))

    # Spacer Storage
    spacerlen = np.zeros(
        len(variantlens),
        dtype=np.int64)

    # Computation Loop
    for idx,vl in enumerate(variantlens):

        # Compute Spacer Length
        spacerlen[idx] += oligolimit - round(vl)

        # Show Update
        liner.send(
            ' Variant {:>{}}: Allows {:,} Base Pair Spacer'.format(
                str(indf.index[idx]),
                plen,
                spacerlen[idx]))

    # Show Time Elapsed
    liner.send('|* Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Results
    return (spacerlen,
        variantlens)

def get_grouped_spacerlen(
    spacerlen,
    liner):
    '''
    Group all spacer lengths by their index
    of occurence. Internal use only.

    :: spacerlen
       type - np.array
       desc - an ordered array of all spacer
              lengths for each variant
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0   = tt.time()
    plen = ut.get_printlen(
        value=len(spacerlen))
    spacergroup = cx.defaultdict(cx.deque)

    # Computation Loop
    for idx,sl in enumerate(spacerlen):

        # Group Spacer
        spacergroup[sl].append(idx)

        # Show Update
        liner.send(
            ' Spacer {:{},d}: Grouped w/ {:,} Other Spacers'.format(
                idx+1,
                plen,
                len(spacergroup[sl])))

    # Show Time Elapsed
    liner.send('|* Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Result
    return spacergroup

# Engine Objective and Helper Functions

def show_update(
    idx,
    plen,
    element,
    motif,
    optstatus,
    optstate,
    inittime,
    terminal,
    liner):
    '''
    Display the current progress in motif
    generation. Internal use only.

    :: element
       type - string
       desc - motif element name, e.g.
              'Motif' or 'Spacer'
    :: motif
       type - string
       desc - a partially explored motif
              sequence path
    :: optstatus
       type - integer
       desc - motif feasibility status
    :: optstate
       type - integer
       desc - feasibility failure state marker
    :: inittime
       type - tt.time
       desc - initial time stamp
    :: terminal
       type - boolean
       desc - if True will terminate update to newline
              otherwise, rewrite previous update
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    if len(motif) >= 30:
        design = motif[:14] + '..' + motif[-14:]
    else:
        design = motif

    liner.send(' Candidate {:{},d}: {} {} is {}{}'.format(
        idx,
        plen,
        element,
        design,
        ['Rejected', 'Provisionally Accepted', 'Accepted'][optstatus],
        ['',
        ' due to Excluded Motif',
        ' due to Edge Effect'][optstate]))

    if terminal:
        liner.send('|* Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - inittime))

def is_exmotif_feasible(
    motif,
    exmotifs,
    exmotifindex):
    '''
    Determine if motif devoid of exmotifs.
    Internal use only.

    :: motif
       type - string
       desc - a partially explored motif
              sequence path
    :: exmotifs
       type - set / None
       desc - set of all excluded motifs
    :: exmotifindex
       type - set
       desc - set of constraint embedded
              exmotif ending indices
    '''

    return ut.is_local_exmotif_feasible(
        seq=motif,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex)

def is_edge_feasible(
    motif,
    motiflen,
    lcseq,
    rcseq,
    edgeeffectlength,
    prefixforbidden,
    suffixforbidden):
    '''
    Determine if motif prefix and suffix
    is forbidden. Internal use only.

    :: motif
       type - string
       desc - a paritally explored motif
              sequence path
    :: motiflen
       type - integer
       desc - full motif sequence length
    :: lcseq
       type - string / None
       desc - left context sequence
    :: rcseq
       type - string / None
       desc - right context sequence
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    :: prefixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    '''

    return ut.is_local_edge_feasible(
        seq=motif,
        seqlen=motiflen,
        lcseq=lcseq,
        rcseq=rcseq,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden)

def motif_objectives(
    motif,
    motiflen,
    exmotifs,
    exmotifindex,
    lcseq,
    rcseq,
    edgeeffectlength,
    prefixforbidden,
    suffixforbidden,
    inittime,
    stats,
    idx,
    plen,
    element,
    liner):
    '''
    Determine if a motif satisfies all
    local objectives. Internal use only.

    :: motif
       type - string
       desc - a paritally explored motif
              sequence path
    :: motiflen
       type - integer
       desc - full motif sequence length
    :: exmotifs
       type - set / None
       desc - set of all excluded motifs
    :: exmotifindex
       type - set / None
       desc - set of constraint embedded
              exmotif ending indices
    :: lcseq
       type - string / None
       desc - left context sequence
    :: rcseq
       type - string / None
       desc - right context sequence
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    :: prefixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    :: inittime
       type - tt.time
       desc - initial time stamp
    :: stats
       type - dict
       desc - primer design stats storage
    :: idx
       type - integer
       desc - context assignment index for
              motif being designed
    :: plen
       type - integer
       plen - target index printing length
    :: element
       type - string
       desc - motif element name, e.g.
              'Motif' or 'Spacer'
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Objective 1: Motif Embedding
    obj1, exmotif = is_exmotif_feasible(
        motif=motif,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex)

    # Objective 1 Failed
    if not obj1:

        # Show Update
        show_update(
            idx=idx,
            plen=plen,
            element=element,
            motif=motif,
            optstatus=0,
            optstate=1,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['exmotif_fail'] += 1
        stats['vars']['exmotif_counter'][exmotif] += 1

        # Return Traceback
        return False, max(
            0,
            len(motif)-1)

    # Objective 2: Edge Feasibility (Edge-Effects)
    obj2, dxmotifs, traceloc = is_edge_feasible(
        motif=motif,
        motiflen=motiflen,
        lcseq=lcseq,
        rcseq=rcseq,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden)

    # Objective 2 Failed
    if not obj2:

        # Show Update
        show_update(
            idx=idx,
            plen=plen,
            element=element,
            motif=motif,
            optstatus=0,
            optstate=2,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['edge_fail'] += len(dxmotifs)
        stats['vars']['exmotif_counter'].update(dxmotifs)

        # Return Traceback
        return False, traceloc

    # Show Update
    show_update(
        idx=idx,
        plen=plen,
        element=element,
        motif=motif,
        optstatus=1,
        optstate=0,
        inittime=inittime,
        terminal=False,
        liner=liner)

    # All Objectives OK!
    return True

def extra_assign_motif(
    motif,
    contextarray,
    leftselector,
    rightselector,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    storage,
    stats,
    plen,
    element,
    element_key,
    liner):
    '''
    Reassign motifs to additional contexts where
    edge effects are absent. Internal use only.

    :: motif
       type - string
       desc - a fully explored motif
              sequence path
    :: contextarray
       type - np.array
       desc - context assignment array
    :: leftselector
       type - lambda
       desc - selector for the left
              sequence context
    :: rightselector
       type - lambda
       desc - selector for the right
              sequence context
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    :: prefixdict
       type - dict
       desc - dictionary of all forbidden
              motif prefixes for left
              context sequences
    :: suffixdict
       type - dict
       desc - dictionary of all forbidden
              motif suffixes for right
              context sequences
    :: storage
       type - list
       desc - list of designed motifs for
              each indexed context
    :: stats
       type - dict
       desc - primer design stats storage
    :: idx
       type - integer
       desc - context assignment index for
              motif being designed
    :: element
       type - string
       desc - motif element name, e.g.
              'Motif' or 'Spacer'
    :: element_key
       type - string
       desc - stats key for element storage
              count
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    i = 0
    t = len(contextarray)

    # Loop through contexts for assignment
    while i < t:

        # Fetch Context
        aidx = contextarray.popleft()

        # Fetch Context Sequences
        lcseq =  leftselector(aidx)
        rcseq = rightselector(aidx)

        # Define Forbidden Prefix and Suffix
        prefixforbidden = prefixdict[lcseq] if not lcseq is None else None
        suffixforbidden = suffixdict[rcseq] if not rcseq is None else None

        # Compute Edge Feasibility
        obj, dxmotifs, _ = ut.is_local_edge_feasible(
            seq=motif,
            seqlen=len(motif),
            lcseq=lcseq,
            rcseq=rcseq,
            edgeeffectlength=edgeeffectlength,
            prefixforbidden=prefixforbidden,
            suffixforbidden=suffixforbidden)

        # Objective Met
        if obj:

            # Record Designed Motif
            storage[aidx] = motif
            stats['vars'][element_key] += 1

            # Show Update
            show_update(
                idx=aidx+1,
                plen=plen,
                element=element,
                motif=motif,
                optstatus=2,
                optstate=0,
                inittime=None,
                terminal=False,
                liner=liner)

        # Objective Failed
        else:

            # Record Failure Stats
            stats['vars']['edge_fail'] += len(dxmotifs)
            stats['vars']['exmotif_counter'].update(dxmotifs)

            # Try Again Later
            contextarray.append(aidx)

        # Update Iteration
        i += 1

def motif_engine(
    motifseq,
    homology,
    optrequired,
    leftcontext,
    rightcontext,
    exmotifs,
    exmotifindex,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    targetcount,
    stats,
    liner):
    '''
    Compute edge-effect free constrained
    motifs within given contexts.

    :: motifseq
       type - string
       desc - degenerate motif design sequence
              constraint
    :: homology
       type - integer
       desc - maximum allowed internal repeat
              length for default traceback
    :: optrequired
       type - bool
       desc - if True, then constraint is degenerate
              so optimization is required
    :: leftcontext
       type - list / None
       desc - list of sequence to the
              left of motifs,
              None otherwise
    :: rightcontext
       type - list / None
       desc - list of sequences to the
              right of motifs,
              None otherwise
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed motifs,
              None otherwise
    :: exmotifindex
       type - set / None
       desc - location of excluded motifs
              in motif design constraint,
              None otherwise
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixdict
       type - dict / None
       desc - dictionary of forbidden motif
              prefix sequences
    :: suffixdict
       type - dict / None
       desc - dictionary of forbidden motif
              suffix sequences
    :: targetcount
       type - integer
       desc - required number of motifs
              to be designed
    :: stats
       type - dict
       desc - motif design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0     = tt.time()          # Start Timer
    motifs = [None]*targetcount # Store Motifs
    plen   = ut.get_printlen(   # Target Print Length
        value=targetcount)

    # Context Setup
    contextarray   = cx.deque(range(targetcount))  # Context Array
    (_,
    leftselector)  = ut.get_context_type_selector( # Left  Context Selector
        context=leftcontext)
    (_,
    rightselector) = ut.get_context_type_selector( # Right Context Selector
        context=rightcontext)

    # Optimize exmotifs
    if not exmotifs is None:
        exmotifs = ut.get_grouped_sequences(
            sequences=exmotifs)

    # Motif Design Required
    if optrequired:

        # Define Maker Instance
        maker = nr.base.maker.NRPMaker(
            part_type='DNA',
            seed=None)

        # Core Design Loop
        while contextarray:

            # Fetch Context Index
            idx = contextarray.popleft()

            # Fetch Context Sequences
            lcseq =  leftselector(idx)
            rcseq = rightselector(idx)

            # Define Forbidden Prefix and Suffix
            prefixforbidden = prefixdict[lcseq] if not lcseq is None else None
            suffixforbidden = suffixdict[rcseq] if not rcseq is None else None

            # Define Objective Function
            objectivefunction = lambda motif: motif_objectives(
                motif=motif,
                motiflen=len(motifseq),
                exmotifs=exmotifs,
                exmotifindex=exmotifindex,
                lcseq=lcseq,
                rcseq=rcseq,
                edgeeffectlength=edgeeffectlength,
                prefixforbidden=prefixforbidden,
                suffixforbidden=suffixforbidden,
                inittime=t0,
                stats=stats,
                idx=idx+1,
                plen=plen,
                element='Motif',
                liner=liner)

            # Design Motif via Maker
            motif = maker.nrp_maker(
                homology=min(len(motifseq), homology),
                seq_constr=motifseq,
                struct_constr='.'*len(motifseq),
                target_size=1,
                background=None,
                struct_type=None,
                synth_opt=False,
                local_model_fn=objectivefunction,
                jump_count=100,
                fail_count=100,
                output_file=None,
                verbose=False,
                abortion=True,
                allow_internal_repeat=True,
                check_constraints=False)

            # Did we succeed? No ..
            if len(motif) == 0:

                # Terminate Design Loop
                break # RIP .. We failed!

            # A motif was designed!
            else:

                # Extract Motif
                motif = motif[0]

                # Record Designed Motif
                motifs[idx] = motif
                stats['vars']['motif_count'] += 1

                # Show Update
                show_update(
                    idx=idx+1,
                    plen=plen,
                    element='Motif',
                    motif=motifs[idx],
                    optstatus=2,
                    optstate=0,
                    inittime=None,
                    terminal=False,
                    liner=liner)

                extra_assign_motif(
                    motif=motif,
                    contextarray=contextarray,
                    leftselector=leftselector,
                    rightselector=rightselector,
                    edgeeffectlength=edgeeffectlength,
                    prefixdict=prefixdict,
                    suffixdict=suffixdict,
                    storage=motifs,
                    stats=stats,
                    plen=plen,
                    element='Motif',
                    element_key='motif_count',
                    liner=liner)

    # Constant Motif
    else:

        # Constant Solution
        motifs = motifseq
        stats['vars']['motif_count'] = targetcount

    # Check Status and Return Solution
    if not optrequired or \
       stats['vars']['motif_count'] == targetcount:

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

        # Determine Last Known Motif
        if optrequired:
            lastmotif = motifs[-1]
        else:
            lastmotif = motifseq

        # Final Update
        show_update(
            idx=targetcount,
            plen=plen,
            element='Motif',
            motif=lastmotif,
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Return Results
        return (motifs, stats)

    # Design Unsuccessful
    else:

        # This was a miscarriage
        stats['status'] = False
        stats['basis']  = 'unsolved'

        # Final Update
        liner.send('|* Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

        # Return Results
        return (None, stats)

def spacer_engine(
    spacergroup,
    leftcontext,
    rightcontext,
    exmotifs,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    targetcount,
    stats,
    liner):
    '''
    Compute constrained degenerate
    spacer contextually.

    :: spacergroup
       type - cx.defaultdict
       desc - grouped spacer lengths by their index
              of occurence
    :: leftcontext
       type - list / None
       desc - list of sequence to the
              left of spacers,
              None otherwise
    :: rightcontext
       type - list / None
       desc - list of sequences to the
              right of spacers,
              None otherwise
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed spacers,
              None otherwise
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixdict
       type - dict / None
       desc - dictionary of forbidden motif
              prefix sequences
    :: suffixdict
       type - dict / None
       desc - dictionary of forbidden motif
              suffix sequences
    :: targetcount
       type - integer
       desc - required number of spacers
              to be designed
    :: stats
       type - dict
       desc - motif design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0      = tt.time()          # Start Timer
    spacers = [None]*targetcount # Store Spacers
    plen    = ut.get_printlen(   # Target Print Length
        value=targetcount)
    abort   = False              # Abortion Flag

    # Context Setup
    contextarray   = cx.deque(range(targetcount))  # Context Array
    (_,
    leftselector)  = ut.get_context_type_selector( # Left  Context Selector
        context=leftcontext)
    (_,
    rightselector) = ut.get_context_type_selector( # Right Context Selector
        context=rightcontext)

    # Optimize exmotifs
    if not exmotifs is None:
        exmotifs = ut.get_grouped_sequences(
            sequences=exmotifs)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Core Design Outer Loop
    while spacergroup:

        # Fetch Context Array
        (spacerlen,
        contextarray) = spacergroup.popitem()

        # Core Design Inner Loop
        while contextarray:

            # Case 1: Zero Length Spacer
            if spacerlen == 0:

                # Assign Gap to Entire Group
                while contextarray:

                    # Fetch Context Index
                    idx = contextarray.popleft()

                    # Assign Gap
                    spacers[idx] = '-'
                    stats['vars']['spacer_count'] += 1

                    # Show Update
                    show_update(
                        idx=idx+1,
                        plen=plen,
                        element='Spacer',
                        motif='*GAP*',
                        optstatus=2,
                        optstate=0,
                        inittime=None,
                        terminal=False,
                        liner=liner)

                # Zero Spacers Completed
                break # Jump to Next Group

            # Case 2: Non-Zero Length Spacer
            #         (this requires computation)
            else:

                # Fetch Context Index
                idx = contextarray.popleft()

                # Fetch Context Sequences
                lcseq =  leftselector(idx)
                rcseq = rightselector(idx)

                # Define Forbidden Prefix and Suffix
                prefixforbidden = prefixdict[lcseq] if not lcseq is None else None
                suffixforbidden = suffixdict[rcseq] if not rcseq is None else None

                # Define Objective Function
                objectivefunction = lambda spacer: motif_objectives(
                    motif=spacer,
                    motiflen=spacerlen,
                    exmotifs=exmotifs,
                    exmotifindex=None,
                    lcseq=lcseq,
                    rcseq=rcseq,
                    edgeeffectlength=edgeeffectlength,
                    prefixforbidden=prefixforbidden,
                    suffixforbidden=suffixforbidden,
                    inittime=t0,
                    stats=stats,
                    idx=idx+1,
                    plen=plen,
                    element='Spacer',
                    liner=liner)

                # Design Spacer via Maker
                spacer = maker.nrp_maker(
                    homology=min(spacerlen, 6),
                    seq_constr='N'*spacerlen,
                    struct_constr='.'*spacerlen,
                    target_size=1,
                    background=None,
                    struct_type=None,
                    synth_opt=False,
                    local_model_fn=objectivefunction,
                    jump_count=100,
                    fail_count=100,
                    output_file=None,
                    verbose=False,
                    abortion=True,
                    allow_internal_repeat=True,
                    check_constraints=False)

                # Did we succeed? No ..
                if len(spacer) == 0:

                    # Terminate Design Loop
                    abort = True
                    break # RIP .. We failed!

                # A spacer was designed!
                else:

                    # Extract Spacer
                    spacer = spacer[0]

                    # Record Designed Motif
                    spacers[idx] = spacer
                    stats['vars']['spacer_count'] += 1

                    # Show Update
                    show_update(
                        idx=idx+1,
                        plen=plen,
                        element='Spacer',
                        motif=spacers[idx],
                        optstatus=2,
                        optstate=0,
                        inittime=None,
                        terminal=False,
                        liner=liner)

                    extra_assign_motif(
                        motif=spacer,
                        contextarray=contextarray,
                        leftselector=leftselector,
                        rightselector=rightselector,
                        edgeeffectlength=edgeeffectlength,
                        prefixdict=prefixdict,
                        suffixdict=suffixdict,
                        storage=spacers,
                        stats=stats,
                        plen=plen,
                        element='Spacer',
                        element_key='spacer_count',
                        liner=liner)

        # Continue Outer Loop?
        if abort:
            break

    # Check Status and Return Solution
    if not abort and \
       stats['vars']['spacer_count'] == targetcount:

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

         # Determine Last Known Motif
        if spacers[-1] != '-':
            lastmotif = spacers[-1]
        else:
            lastmotif = '*GAP*'

        # Final Update
        show_update(
            idx=targetcount,
            plen=plen,
            element='Spacer',
            motif=lastmotif,
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Return Results
        return (spacers, stats)

    # Design Unsuccessful
    else:

        # This was a miscarriage
        stats['status'] = False
        stats['basis']  = 'unsolved'

        # Final Update
        liner.send('|* Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

        # Return Results
        return (None, stats)
