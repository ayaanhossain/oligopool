import os
import time as tt

from .base import utils as ut
from .base import core_inspect as ci


def inspect(
    target:str,
    kind:str='auto',
    verbose:bool=True) -> dict:
    '''
    Inspect non-CSV artifacts produced by Oligopool Calculator (background DBs, index files, pack files).

    Required Parameters:
        - `target` (`str`): Path to an artifact (background directory, `.oligopool.index`, or `.oligopool.pack`).

    Optional Parameters:
        - `kind` (`str`): One of `background`, `index`, `pack`, or `auto` (default: `auto`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A dictionary of stats containing a compact summary of the inspected artifact.

    Notes:
        - `inspect` is read-only: it does not repair artifacts or deserialize unsafe objects (no unpickling).
        - Background directories are validated before opening to avoid accidentally creating a new DB.
    '''

    # Argument Aliasing
    artifactpath = target
    kind         = kind
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Inspect Verbiage Print
    liner.send('\n[Oligopool Calculator: QC Mode - Inspect]\n')

    # Start Timer
    t0 = tt.time()

    # Init Stats
    stats = {
        'status'   : False,
        'basis'    : 'infeasible',
        'step'     : 1,
        'step_name': 'inspecting-artifact',
        'vars'     : None,
        'warns'    : {},
    }

    # Validate target
    if not isinstance(artifactpath, str) or not artifactpath:
        liner.send(
            '\n Target: {} [INVALID]\n'.format(artifactpath))
        liner.close()
        return ut.stamp_stats(
            stats=stats, module='inspect',
            input_rows=0, output_rows=0)

    # Normalize and validate kind
    if kind is None:
        kind_norm = 'auto'
    elif not isinstance(kind, str):
        kind_norm = None
    else:
        kind_norm = kind.strip().lower()
    if kind_norm in ('', 'none'):
        kind_norm = 'auto'
    if kind_norm not in ('auto', 'background', 'index', 'pack'):
        liner.send(
            '\n Target: {}\n'.format(artifactpath))
        liner.send(
            '   Kind: Invalid "{}"; use background, index, or pack.\n'.format(
                kind))
        liner.close()
        return ut.stamp_stats(
            stats=stats, module='inspect',
            input_rows=0, output_rows=0)

    # Detect kind conservatively (avoid accidental DB creation)
    detected = None if kind_norm == 'auto' else kind_norm
    if detected is None:
        if artifactpath.endswith('.oligopool.index'):
            detected = 'index'
        elif artifactpath.endswith('.oligopool.pack'):
            detected = 'pack'
        elif ci.is_vectordb_dir(artifactpath) or \
             artifactpath.endswith('.oligopool.background'):
            detected = 'background'

    if detected is None:
        liner.send(
            '\n Target: {}\n'.format(artifactpath))
        liner.send(
            '   Kind: Could not infer; use kind=background, index, or pack.\n')
        liner.close()
        return ut.stamp_stats(
            stats=stats, module='inspect',
            input_rows=0, output_rows=0)

    # Human-readable kind labels
    _KIND_LABEL = {
        'background': 'Background Database',
        'index'     : 'Index File',
        'pack'      : 'Pack File',
    }

    # Print header
    inferred = kind_norm == 'auto'
    liner.send(
        '\n Target: {}\n'.format(artifactpath))
    liner.send(
        '   Kind: Target {} as {}\n'.format(
            'Inferred' if inferred else 'Evaluated',
            _KIND_LABEL[detected]))

    # Inspect artifact
    ok    = False
    warns = {}
    meta  = {}

    if detected == 'background':
        norm_dir = ci.normalize_background_dir(artifactpath)
        if not os.path.isdir(norm_dir):
            liner.send(
                '\n [ERROR] Directory does not exist: {}\n'.format(
                    norm_dir))
        elif not ci.is_vectordb_dir(artifactpath):
            liner.send(
                '\n [ERROR] Not a valid background DB: {}\n'.format(
                    norm_dir))
        else:
            ok, warns, meta = ci.inspect_background(artifactpath)

    elif detected == 'index':
        resolved = ut.get_adjusted_path(artifactpath, '.oligopool.index')
        if not os.path.isfile(resolved):
            liner.send(
                '\n [ERROR] File does not exist: {}\n'.format(
                    resolved))
        else:
            ok, warns, meta = ci.inspect_zip(
                path=resolved,
                load_dict_names=('meta.map',))

    elif detected == 'pack':
        resolved = ut.get_adjusted_path(artifactpath, '.oligopool.pack')
        if not os.path.isfile(resolved):
            liner.send(
                '\n [ERROR] File does not exist: {}\n'.format(
                    resolved))
        else:
            ok, warns, meta = ci.inspect_zip(
                path=resolved,
                load_dict_names=('packing.stat',))

    # Inspection Summary
    liner.send('\n[Inspection Summary]\n')

    if ok:
        if detected == 'background':
            K   = meta.get('K', 0)
            LEN = meta.get('LEN', 0)
            liner.send(
                ' Background Repeat Length: {:,}\n'.format(K - 1))
            liner.send(
                ' Background   Fill Count : {:,}\n'.format(LEN))

        elif detected == 'index':
            entries     = meta.get('entries', [])
            total_bytes = sum(e.get('bytes', 0) for e in entries)
            mm = meta.get('meta.map', {}) if isinstance(
                meta.get('meta.map'), dict) else {}
            liner.send(
                '  Archive Entries: {:,}\n'.format(
                    len(entries)))
            liner.send(
                '    Total Size   : {:,} Byte(s)\n'.format(
                    total_bytes))
            if 'barcodename' in mm:
                liner.send(
                    '  Barcode Column : {}\n'.format(
                        mm['barcodename']))
            if 'variantcount' in mm:
                liner.send(
                    '  Barcode Count  : {:,}\n'.format(
                        mm['variantcount']))
            if 'barcodelen' in mm:
                liner.send(
                    '  Barcode Length : {:,} bp\n'.format(
                        mm['barcodelen']))
            if mm.get('barcodeprefix'):
                liner.send(
                    '   Prefix Anchor : {} ({:,} bp)\n'.format(
                        mm['barcodeprefix'], len(mm['barcodeprefix'])))
                if mm.get('barcodepregap'):
                    liner.send(
                        '   Prefix Gap    : {:,} bp\n'.format(
                            mm['barcodepregap']))
            if mm.get('barcodesuffix'):
                liner.send(
                    '   Suffix Anchor : {} ({:,} bp)\n'.format(
                        mm['barcodesuffix'], len(mm['barcodesuffix'])))
                if mm.get('barcodepostgap'):
                    liner.send(
                        '   Suffix Gap    : {:,} bp\n'.format(
                            mm['barcodepostgap']))
            liner.send(
                '     Association : {}\n'.format(
                    'Yes' if mm.get('association') else 'No'))
            if mm.get('associateprefix'):
                liner.send(
                    '   Assoc. Prefix : {} ({:,} bp)\n'.format(
                        mm['associateprefix'], len(mm['associateprefix'])))
                if mm.get('associatepregap'):
                    liner.send(
                        '   Assoc. PreGap : {:,} bp\n'.format(
                            mm['associatepregap']))
            if mm.get('associatesuffix'):
                liner.send(
                    '   Assoc. Suffix : {} ({:,} bp)\n'.format(
                        mm['associatesuffix'], len(mm['associatesuffix'])))
                if mm.get('associatepostgap'):
                    liner.send(
                        '   Assoc. PostGap: {:,} bp\n'.format(
                            mm['associatepostgap']))

        elif detected == 'pack':
            entries     = meta.get('entries', [])
            total_bytes = sum(e.get('bytes', 0) for e in entries)
            pstat = meta.get('packing.stat', {}) if isinstance(
                meta.get('packing.stat'), dict) else {}
            _PACK_TYPES = {0: 'concatenate', 1: 'merge'}
            plen = ut.get_printlen(
                value=max(pstat.get('scanned_reads', 0),
                          pstat.get('survived_reads', 0),
                          pstat.get('packed_reads', 0)))
            liner.send(
                '  Archive Entries: {:,}\n'.format(
                    len(entries)))
            liner.send(
                '    Total Size   : {:,} Byte(s)\n'.format(
                    total_bytes))
            if 'pack_type' in pstat:
                liner.send(
                    '     Pack Type   : {}\n'.format(
                        _PACK_TYPES.get(pstat['pack_type'], pstat['pack_type'])))
            if 'pack_size' in pstat:
                liner.send(
                    '     Pack Size   : {:.1f}M Reads/Pack\n'.format(
                        pstat['pack_size']))
            if 'pack_count' in pstat:
                liner.send(
                    '     Pack Count  : {:,}\n'.format(
                        pstat['pack_count']))
            if 'scanned_reads' in pstat:
                liner.send(
                    '  Scanned Reads  : {:{},d}\n'.format(
                        pstat['scanned_reads'], plen))
            if 'survived_reads' in pstat:
                liner.send(
                    ' Survived Reads  : {:{},d}\n'.format(
                        pstat['survived_reads'], plen))
            if 'packed_reads' in pstat:
                liner.send(
                    '   Packed Reads  : {:{},d}\n'.format(
                        pstat['packed_reads'], plen))

    # Show Verdict
    corrupt = bool(warns.get('missing_entries') or warns.get('failed_load'))
    if ok and not corrupt:
        verdict = 'Valid'
    elif ok and corrupt:
        verdict = 'Corrupted'
    else:
        verdict = 'Invalid'
    liner.send(
        ' Verdict: {} {}\n'.format(
            _KIND_LABEL.get(detected, 'Artifact'), verdict))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Build Return Stats
    stats['status'] = bool(ok) and not corrupt
    if ok and not corrupt:
        stats['basis'] = 'solved'
    elif ok and corrupt:
        stats['basis'] = 'corrupted'
    else:
        stats['basis'] = 'infeasible'
    stats['vars']   = {
        'kind'   : detected,
        'target' : artifactpath,
        'verdict': verdict,
        'meta'   : meta,
    }
    stats['warns'] = warns

    # Close Liner
    liner.close()

    # Return Statistics
    return ut.stamp_stats(
        stats=stats, module='inspect',
        input_rows=0, output_rows=0)
