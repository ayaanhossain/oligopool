import os
import time as tt

import collections as cx
import zipfile     as zf
import random      as rn
import atexit      as ae

import ctypes                       as ct
import multiprocessing              as mp
import multiprocessing.sharedctypes as mpsct

import numpy  as np
import pandas as pd
import edlib  as ed

import utils     as ut
import valparse  as vp
import scry      as sy
import phiX      as px

def count_engine(
    prepfile,
    packfile,
    packqueue,
    countqueue,
    amperrors,
    barerrors,
    coreid,
    ncores,
    analyzedreads,
    diagnosticreads,
    phiXreads,
    mononukereads,
    dinukereads,
    trinukereads,
    incalcreads,
    liner):
    '''
    Count read packs stored in packfile and
    enqueue the final read count matrix in to
    countqueue. Internal use only.

    :: prepfile
       type - string
       desc - path to compressed zipfile
              storing prepared structures
              and data models for amplicon
              and sample data
    :: packfile
       type - string
       desc - path to compressed zipfile
              storing read packs
    :: packqueue
       type - SimpleQueue
       desc - queue storing path to read pack
              files stored in packfile
    :: countqueue
       type - SimpleQueue
       desc - queue storing count matrices
              aggregating one or more read
              packs processed by each core
    :: amperrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in amplicons / controls
              before discarding reads from
              counting
    :: barerrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in barcodes before
              discarding reads from counting
    :: coreid
       type - integer
       desc - current core integer id
    :: ncores
       type - integer
       desc - total number of counters
              concurrently initiated
    :: analyzedreads
       type - c_longlong
       desc - total number of existing read
              pairs in packfile
    :: diagnosticreads
       type - c_longlong
       desc - total number of read pairs
              of sample origin
    :: phiXreads
       type - c_longlong
       desc - total number of read pairs
              discarded of PhiX origin
    :: phiXreads
       type - c_longlong
       desc - total number of read pairs
              discarded due to dinucleotide
              read composition
    :: mononukereads
       type - c_longlong
       desc - total number of read pairs
              discarded and containing
              homopolymer runs
    :: dinukereads
       type - c_longlong
       desc - total number of read pairs
              discarded and containing
              dinucleotide runs
    :: trinukereads
       type - c_longlong
       desc - total number of read pairs
              discarded and containing
              dinucleotide runs
    :: incalcreads
       type - c_longlong
       desc - total number of read pairs
              discarded due to one or more
              missing elements, and of
              unknown origin (bad products)
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Open prepfile
    prepfile = zf.ZipFile(
        file=prepfile)

    # Load Models
    wrmodel = ut.loadmodel(
        archive=prepfile,
        mfile='wr.model')
    prmodel = ut.loadmodel(
        archive=prepfile,
        mfile='pr.model')
    acmodel = ut.loadmodel(
        archive=prepfile,
        mfile='ac.model')

    # Load Maps
    metamap = ut.loaddict(
        archive=prepfile,
        dfile='meta.map')
    sammap = ut.loaddict(
        archive=prepfile,
        dfile='sam.map')

    # Close prepfile
    prepfile.close()

    # Adjust Errors
    amperrors = min(
        metamap['actval'],
        amperrors)
    barerrors = min(
        min(
            metamap['prtval'],
            metamap['wrtval']),
        barerrors)

    # Define PhiX Spectrum
    phiXkval = 30
    phiXspec = px.get_phiX_spectrum(
        k=phiXkval)

    # Open packfile
    packfile = zf.ZipFile(
        file=packfile)

    # Load packing.stat
    packstat = ut.loaddict(
        archive=packfile,
        dfile='packing.stat')

    # Extract Sequencing Type
    seqtype = packstat['seqtype']

    # Book-keeping Variables
    c_analyzedreads   = 0
    c_diagnosticreads = 0
    c_phiXreads       = 0
    c_mononukereads   = 0
    c_dinukereads     = 0
    c_trinukereads    = 0
    c_incalcreads     = 0

    # Compute Printing Lengths
    clen = len(str(ncores))
    plen = ut.get_printlen(
        value=packstat['packsize'])
    slen = plen = ut.get_printlen(
        value=packstat['survived'])

    # Read Count Matrix (Storage)
    countmatrix = np.zeros((
        metamap['samsize'],
        len(metamap['cfields'])),
        dtype=np.float64)

    # Initial Update
    if packqueue.empty():
        liner.send(
            ' Core {:{},d}: No Packs to Analyze ... Exiting'.format(
                coreid,
                clen))

    # Pack Counting Loop
    while not packqueue.empty():

        # Fetch Pack Path / Token
        cpath = packqueue.get()

        # Exhaustion Token
        if cpath is None:
            break

        # Fetch Read Pack
        cpack = ut.loadpack(
            archive=packfile,
            pfile=cpath)

        # Book-keeping Variables
        exoread = None
        exofreq = None

        verbagereach  = 0
        verbagetarget = rn.randint(
            *map(round, (len(cpack) * 0.080,
                         len(cpack) * 0.120)))

        readcount = 0

        packname  = cpath.rstrip('.pack')

        # Start Timer
        t0 = tt.time()

        # Read Counting Loop
        while True:

            # Exoneration Block
            if not exoread is None:

                # Exoneration Flag
                exonerated = False

                # PhiX Match?
                for kmer in ut.stream_spectrum(
                    seq=exoread,
                    k=phiXkval):

                    if kmer in phiXspec:
                        c_phiXreads += exofreq
                        exonerated   = True
                        break

                # Nucleotide Composition
                acount = exoread.count('A')
                tcount = exoread.count('T')
                gcount = exoread.count('G')
                ccount = exoread.count('C')

                # Dinucleotide Match?
                if not exonerated:
                    # Composition Largely Dinucleotide?
                    dinukethresh = 0.75 * len(exoread)
                    if (ccount + tcount >= dinukethresh) or \
                       (acount + ccount >= dinukethresh) or \
                       (ccount + gcount >= dinukethresh) or \
                       (acount + gcount >= dinukethresh) or \
                       (acount + tcount >= dinukethresh) or \
                       (tcount + gcount >= dinukethresh):
                        c_dinukereads += exofreq
                        exonerated = True
                        # print(exoread)

                # Mononucleotide Match?
                if not exonerated:
                    # Composition Large Mononucleotide?
                    mononukethresh = 0.50 * len(exoread)
                    if ccount >= mononukethresh or \
                       acount >= mononukethresh or \
                       gcount >= mononukethresh or \
                       tcount >= mononukethresh:
                        c_mononukereads += exofreq
                        exonerated = True

                # Trinucleotide Match?
                if not exonerated:
                    # Composition Largely Dinucleotide?
                    trinukethresh = 1.00 * len(exoread)
                    if (ccount + acount + tcount >= trinukethresh) or \
                       (ccount + acount + gcount >= trinukethresh) or \
                       (acount + tcount + gcount >= trinukethresh) or \
                       (tcount + ccount + gcount >= trinukethresh):
                        c_dinukereads += exofreq
                        exonerated = True
                        # print(exoread)

                # Uncategorized Discard
                if not exonerated:
                    c_incalcreads += exofreq

                # Clear for Next Exoneration
                exoread = None
                exofreq = None

            # Time to Show Update?
            if verbagereach >= verbagetarget:

                # Show Update
                liner.send(
                    ' Core {:{},d}: Analyzed {:{},d} Reads in {:.2f} sec'.format(
                        coreid,
                        clen,
                        c_analyzedreads,
                        slen,
                        tt.time()-t0))

                # Update Book-keeping
                verbagereach = 0

            # Continue processing pack?
            if not cpack:
                break

            # Fetch Read and Frequency
            # ((ampread,
            # barread),
            # crpfreq) = cpack.pop()

            # Fetch Pack Entry
            cpentry = cpack.pop()

            # Parse Pack Entry
            if seqtype == 2: # Paired-End Sequencing
                ((ampread,
                barread),
                crpfreq) = cpentry
            else:            # Single-End Sequencing
                (fullread,
                crpfreq) = cpentry
                ampread  = fullread
                barread  = ut.get_revcomp(
                    seq=ampread)

            # Update Book-keeping
            c_analyzedreads += crpfreq
            verbagereach    += 1
            readcount       += 1

            # Locate revconst
            aln = ed.align(
                query=metamap['revconst'],
                target=barread,
                mode='HW',
                task='locations',
                k=2)

            # revconst Absent
            if aln['editDistance'] == -1:
                exoread = ampread
                exofreq = crpfreq
                continue

            # Determine PlateReverseBarcode ID
            prend = aln['locations'][-1][+0] + 0
            prstr = aln['locations'][+0][+0] + 0 - \
                    metamap['prlen']
            prid  = prmodel.predict(
                x=barread[prstr:prend],
                mx=barerrors)

            # Plate Reverse Barcode Absent
            if prid is None:
                c_incalcreads += crpfreq
                continue

            # Determine Well Reverse Barcode
            wrstr = aln['locations'][+0][-1] + 1
            wrend = aln['locations'][-1][-1] + 1 + \
                    metamap['wrlen']
            wrid  = wrmodel.predict(
                x=barread[wrstr:wrend],
                mx=barerrors)

            # Well Reverse Barcode Absent
            if wrid is None:
                c_incalcreads += crpfreq
                continue

            # (prid, wrid) Combination Absent
            if not (prid, wrid) in sammap:
                c_incalcreads += crpfreq
                continue

            # Determine Amplicon / Control
            aid = acmodel.predict(
                x=ampread,
                mx=amperrors)

            # Amplicon / Control Absent
            if aid is None:
                exoread = ampread
                exofreq = crpfreq
                continue

            # All Components Valid
            sid = sammap[(prid, wrid)]
            cid = metamap['cfields'][aid]
            countmatrix[sid, cid] += crpfreq
            c_diagnosticreads     += crpfreq

            # # Time to show updates?
            # if verbagereach >= verbagetarget:

            #     # Show Update
            #     liner.send(
            #         ' Core {:{},d}: Counted Pack {} w/ {:{},d} Reads in {:.2f} sec'.format(
            #             coreid,
            #             clen,
            #             packname,
            #             readcount,
            #             plen,
            #             tt.time()-t0))

            #     # Update Book-keeping
            #     verbagereach = 0

        # Show Final Updates
        liner.send(
            ' Core {:{},d}: Counted Pack {} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                coreid,
                clen,
                packname,
                readcount,
                plen,
                tt.time()-t0))

    # Close packfile
    packfile.close()

    # Update countqueue
    countqueue.put(countmatrix)
    countqueue.put(None)

    # Update Read Counting Book-keeping
    analyzedreads.increment(incr=c_analyzedreads)
    diagnosticreads.increment(incr=c_diagnosticreads)
    phiXreads.increment(incr=c_phiXreads)
    mononukereads.increment(incr=c_mononukereads)
    dinukereads.increment(incr=c_dinukereads)
    trinukereads.increment(incr=c_trinukereads)
    incalcreads.increment(incr=c_incalcreads)

def acount(
    indexfile,
    packfile,
    countfile,
    maptype=0,
    barcodeerrors=-1,
    associateerrors=-1,
    ncores=0,
    memlimit=0,
    verbose=True):
    '''
    TBD
    '''

    # Start Liner
    liner = ut.liner_engine(online=verbose)

    # Packing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Associate Count]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full indexfile Validation
    indexfile_valid = vp.get_indexfile_validity(
        indexfile=indexfile,
        indexfile_field='     Index File  ',
        associated=True,
        liner=liner)

    # Adjust indexfile Suffix
    if indexfile_valid:
        indexfile = ut.get_adjusted_path(
            path=indexfile,
            suffix='.oligopool.index')

    # Full packfile Validation
    (packfile_valid,
    packcount) = vp.get_parsed_packfile(
        packfile=packfile,
        packfile_field='      Pack File  ',
        liner=liner)

    # Adjust packfile Suffix
    if packfile_valid:
        packfile = ut.get_adjusted_path(
            path=packfile,
            suffix='.oligopool.pack')

    # Full countfile Validation
    countfile_valid = vp.get_outfile_validity(
        outfile=countfile,
        outfile_suffix='.oligopool.acount.csv',
        outfile_field='     Count File  ',
        liner=liner)

    # Adjust countfile Suffix
    if not countfile is None:
        countfile = ut.get_adjusted_path(
            path=countfile,
            suffix='.oligopool.acount.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full maptype Validation
    maptype_valid = vp.get_categorical_validity(
        category=maptype,
        category_field='   Mapping Type  ',
        category_pre_desc=' ',
        category_post_desc=' Classification',
        category_dict={
            0: 'Terminus Optimized',
            1: 'Spectrum Optimized'},
        liner=liner)

    # Full barcodeerrors Validation
    barcodeerrors_valid = vp.get_errors_validity(
        errors=barcodeerrors,
        errors_field='   Barcode Errors',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Barcode',
        errors_base='B',
        indexfile_valid=indexfile_valid,
        indexfile=indexfile,
        liner=liner)

    # Full associateerrors Validation
    associateerrors_valid = vp.get_errors_validity(
        errors=associateerrors,
        errors_field=' Associate Errors',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Associate',
        errors_base='A',
        indexfile_valid=indexfile_valid,
        indexfile=indexfile,
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='       Num Cores ',
        default=packcount,
        offset=2,
        liner=liner)

    # Full num_core Parsing and Validation
    (memlimit,
    memlimit_valid) = vp.get_parsed_memory_info(
        memlimit=memlimit,
        memlimit_field='       Mem Limit ',
        ncores=ncores,
        ncores_valid=ncores_valid,
        liner=liner)

    # First Pass Validation
    if not all([
        indexfile_valid,
        packfile_valid,
        countfile_valid,
        maptype_valid,
        barcodeerrors_valid,
        associateerrors_valid,
        ncores_valid,
        memlimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    return

    # Schedule countfile deletion
    ctdeletion = ae.register(
        ut.remove_file,
        countfile)

    # Adjust Errors
    amperrors = round(amperrors)
    barerrors = round(barerrors)

    # Define packqueue
    packqueue = mp.SimpleQueue()

    # Load packqueue
    archive = zf.ZipFile(
        file=packfile)
    for cpath in archive.namelist():
        if cpath.endswith('.pack'):
            packqueue.put(cpath)
    archive.close()
    for _ in range(ncores):
        packqueue.put(None)

    # Define countqueue
    countqueue = mp.SimpleQueue()

    # Read Counting Book-keeping
    analyzedreads   = ut.SafeCounter()
    diagnosticreads = ut.SafeCounter()
    phiXreads       = ut.SafeCounter()
    mononukereads   = ut.SafeCounter()
    dinukereads     = ut.SafeCounter()
    trinukereads    = ut.SafeCounter()
    incalcreads     = ut.SafeCounter()

    # Launching Read Counting
    liner.send('\n[Counting Read Packs]\n')

    # Start Timer
    t0 = tt.time()

    # Counter Process Store
    readcounters = []

    # Fire off Single-Core Read Counter
    if ncores == 1:
        count_engine(
            prepfile=prepfile,
            packfile=packfile,
            packqueue=packqueue,
            countqueue=countqueue,
            amperrors=amperrors,
            barerrors=barerrors,
            coreid=0,
            ncores=ncores,
            analyzedreads=analyzedreads,
            diagnosticreads=diagnosticreads,
            phiXreads=phiXreads,
            mononukereads=mononukereads,
            dinukereads=dinukereads,
            trinukereads=trinukereads,
            incalcreads=incalcreads,
            liner=liner)

    # Fire off Multi-Core Read Counters
    else:

        # Queue All Read Counters
        coreid = 0
        while coreid < ncores:

            # Define Counter
            readcounter = mp.Process(
                target=count_engine,
                args=(prepfile,
                    packfile,
                    packqueue,
                    countqueue,
                    amperrors,
                    barerrors,
                    coreid,
                    ncores,
                    analyzedreads,
                    diagnosticreads,
                    phiXreads,
                    mononukereads,
                    dinukereads,
                    trinukereads,
                    incalcreads,
                    liner,))

            # Start Counter
            readcounter.start()

            # Update Book-keeping
            readcounters.append(readcounter)
            coreid += 1

    # Aggregate All Count Matrices
    agdf = count_aggregator(
        prepfile=prepfile,
        countqueue=countqueue,
        countfile=countfile,
        prodcount=ncores)

    # Join All Read Counters
    for readcounter in readcounters:
        readcounter.join()

    # Counting Status
    if diagnosticreads.value() > 0:
        countstatus = (True,  'Successful')
    else:
        countstatus = (False, 'Failed')

    # Read Counting Stats
    liner.send('\n[Counting Stats]\n')

    plen = ut.get_printlen(
        value=analyzedreads.value())

    liner.send(
        '       Counting  Status: {}\n'.format(
            countstatus[1]))
    liner.send(
        '       Analyzed   Reads: {:{},d}\n'.format(
            analyzedreads.value(),
            plen))
    liner.send(
        '           PhiX   Reads: {:{},d} ({:6.2f} %)\n'.format(
            phiXreads.value(),
            plen,
            ut.safediv(
                A=100. * phiXreads.value(),
                B=analyzedreads.value())))
    liner.send(
        ' Mononucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            mononukereads.value(),
            plen,
            ut.safediv(
                A=100. * mononukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '   Dinucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            dinukereads.value(),
            plen,
            ut.safediv(
                A=100. * dinukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '  Trinucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            trinukereads.value(),
            plen,
            ut.safediv(
                A=100. * trinukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '   Incalculable   Reads: {:{},d} ({:6.2f} %)\n'.format(
            incalcreads.value(),
            plen,
            ut.safediv(
                A=100. * incalcreads.value(),
                B=analyzedreads.value())))
    liner.send(
        '     Diagnostic   Reads: {:{},d} ({:6.2f} %)\n'.format(
            diagnosticreads.value(),
            plen,
            ut.safediv(
                A=(100. * diagnosticreads.value()),
                B=analyzedreads.value())))
    liner.send(
        '           Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Unschedule countfile deletion
    if countstatus[0]:
        ae.unregister(ctdeletion)

    # Close Liner
    liner.close()

    # Return Aggregate DataFrame
    return agdf