import os
import difflib

import numbers      as nu
import multiprocess as mp
import zipfile      as zf
import math         as mt
import collections  as cx

import numpy    as np
import pandas   as pd
import psutil   as pu

from . import vectordb as db
from . import utils    as ut

# Type Parameter Registries
TYPE_REGISTRIES = {
    'barcode_type': {
        'aliases': {
            0: 0, 1: 1,
            'terminus': 0, 'terminus optimized': 0, 'term': 0, 't': 0, 'fast': 0,
            'spectrum': 1, 'spectrum optimized': 1, 'spec': 1, 's': 1, 'slow': 1,
        },
        'display': {0: 'Terminus Optimized', 1: 'Spectrum Optimized'},
        'canonical': ['terminus', 'spectrum'],
    },
    'primer_type': {
        'aliases': {
            0: 0, 1: 1,
            'forward': 0, 'fwd': 0, 'f': 0,
            'reverse': 1, 'rev': 1, 'r': 1,
        },
        'display': {0: 'Forward', 1: 'Reverse'},
        'canonical': ['forward', 'reverse'],
    },
    'motif_type': {
        'aliases': {
            0: 0, 1: 1,
            'per-variant': 0, 'pervariant': 0, 'non-constant': 0, 'nonconstant': 0, 'variable': 0, 'var': 0,
            'constant': 1, 'const': 1, 'anchor': 1, 'fixed': 1,
        },
        'display': {0: 'Non-Constant', 1: 'Constant'},
        'canonical': ['per-variant', 'constant', 'non-constant', 'anchor'],
    },
    'read_type': {
        'aliases': {
            0: 0, 1: 1,
            'forward': 0, 'fwd': 0, 'f': 0,
            'reverse': 1, 'rev': 1, 'r': 1,
        },
        'display': {0: 'Forward', 1: 'Reverse'},
        'canonical': ['forward', 'reverse'],
    },
    'pack_type': {
        'aliases': {
            0: 0, 1: 1,
            'concatenate': 0, 'concatenated': 0, 'concat': 0, 'cat': 0, 'joined': 0, 'join': 0,
            'merge': 1, 'merged': 1, 'assemble': 1, 'assembled': 1, 'asm': 1,
        },
        'display': {0: 'Store Concatenated / Joined Reads', 1: 'Store Assembled / Merged Reads'},
        'canonical': ['concatenate', 'merge'],
    },
    'mapping_type': {
        'aliases': {
            0: 0, 1: 1,
            'fast': 0, 'quick': 0, 'near-exact': 0,
            'sensitive': 1, 'sens': 1, 'accurate': 1, 'slow': 1,
        },
        'display': {0: 'Fast / Near-Exact', 1: 'Slow / Sensitive'},
        'canonical': ['fast', 'sensitive'],
    },
    'join_policy': {
        'aliases': {
            0: 0, 1: 1,
            'left': 0, 'l': 0,
            'right': 1, 'r': 1,
        },
        'display': {0: 'Left-Biased', 1: 'Right-Biased'},
        'canonical': ['left', 'right'],
    },
}


def _normalize_field(field):
    '''
    Normalize field padding for display.
    '''

    return field

def _display_path(path):
    '''
    Return a compact basename-style label for paths.
    Internal use only.
    '''

    if not isinstance(path, str):
        return path
    compact = os.path.basename(path.rstrip(os.sep))
    return compact if compact else path

def get_infile_validity(
    infile,
    infile_suffix,
    infile_field,
    liner):
    '''
    Determine if a given infile exists and
    non-empty. Internal use only.

    :: infile
       type - string
       desc - an input file to check for
              existence and emptiness
    :: infile_suffix
       type - string
       desc - required infile suffix
    :: infile_field
       type - string
       desc - infile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    infile_field = _normalize_field(infile_field)

    # Use space separator for numbered list items (field ends with ])
    sep = '' if infile_field.rstrip().endswith(']') else ':'

    # infile exists and non-empty?
    infile_status = ut.get_path_status(
        path=infile,
        suffix=infile_suffix,
        readable=True,
        writable=False,
        creatable=False)

    # infile is invalid
    if   infile_status == 0:
        liner.send(
            '{}{} {} [INPUT TYPE IS INVALID]\n'.format(
                infile_field, sep, _display_path(infile)))

    else:

        infile = ut.get_adjusted_path(
            path=infile,
            suffix=infile_suffix)

        if   infile_status == 1:
            liner.send(
                '{}{} {} [READ PERMISSION DENIED]\n'.format(
                    infile_field, sep, _display_path(infile)))

        elif infile_status == 3:
            liner.send(
                '{}{} {} [FILE IS EMPTY]\n'.format(
                    infile_field, sep, _display_path(infile)))

        elif infile_status == 'X':
            liner.send(
                '{}{} {} [FILE IS SPECIAL]\n'.format(
                    infile_field, sep, _display_path(infile)))

        elif 5 <= infile_status <= 8:
            liner.send(
                '{}{} {} [FILE IS DIRECTORY]\n'.format(
                    infile_field, sep, _display_path(infile)))

        elif infile_status == 9:
            liner.send(
                '{}{} {} [FILE DOES NOT EXIST]\n'.format(
                    infile_field, sep, _display_path(infile)))

    # infile valid
    return infile_status == 4

def get_indir_validity(
    indir,
    indir_suffix,
    indir_field,
    liner):
    '''
    Determine if a given indir exists and
    non-empty. Internal use only.

    :: indir
       type - string
       desc - an input directory to check for
              existence and emptiness
    :: indir_suffix
       type - string
       desc - required indir suffix
    :: indir_field
       type - string
       desc - indir fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    indir_field = _normalize_field(indir_field)

    # indir exists and non-empty?
    indir_status = ut.get_path_status(
        path=indir,
        suffix=indir_suffix,
        readable=True,
        writable=False,
        creatable=False)

    # indir is invalid
    if   indir_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                indir_field, _display_path(indir)))

    else:

        indir = ut.get_adjusted_path(
            path=indir,
            suffix=indir_suffix)

        if   indir_status == 5:
            liner.send(
                '{}: {} [READ PERMISSION DENIED]\n'.format(
                    indir_field, _display_path(indir)))

        elif indir_status == 7:
            liner.send(
                '{}: {} [DIRECTORY IS EMPTY]\n'.format(
                    indir_field, _display_path(indir)))

        elif indir_status == 'X':
            liner.send(
                '{}: {} [DIRECTORY IS SPECIAL]\n'.format(
                    indir_field, _display_path(indir)))

        elif 1 <= indir_status <= 4:
            liner.send(
                '{}: {} [DIRECTORY IS FILE]\n'.format(
                    indir_field, _display_path(indir)))

        elif indir_status == 9:
            liner.send(
                '{}: {} [DIRECTORY DOES NOT EXIST]\n'.format(
                    indir_field, _display_path(indir)))

    # indir valid
    return indir_status == 8

def get_readfile_validity(
    readfile,
    readfile_field,
    paired_readfile,
    liner):
    '''
    Determine if a given readfile exists, is
    non-empty and of FastQ type? If provided,
    check if given readfile is a duplicate of
    paired_readfile? Internal use only.

    :: readfile
       type - string
       desc - path to FastQ file storing
              reads
    :: readfile_field
       type - string
       desc - readfile fieldname used in
              printing
    :: paired_readfile
       type - string / None
       desc - path to paired FastQ file
              storing reads
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # readfile exists and non-empty?
    readfile_exists = get_infile_validity(
        infile=readfile,
        infile_suffix=None,
        infile_field=readfile_field,
        liner=liner)

    # readfile is of FastQ type?
    readfile_is_fastq = False
    if readfile_exists:

        # FastQ file read attempt
        try:
            next(ut.stream_fastq_engine(
                filepath=readfile))

        # Read unsuccesful
        except:
            liner.send(
                '{}: {} [INVALID FASTQ FILE]\n'.format(
                    readfile_field, _display_path(readfile)))

        # Read successful
        else:
            readfile_is_fastq = True

    # Pairs duplicate?
    readfile_duplicate = False
    if readfile_is_fastq:

        # Pair to compare?
        if not paired_readfile is None:
            if os.path.samefile(
                readfile,
                paired_readfile):
                liner.send(
                    '{}: {} [DUPLICATE OF R1 FILE]\n'.format(
                        readfile_field, _display_path(readfile)))
                readfile_duplicate = True

        # Pair comparison successful
        # or unnecessary
        if not readfile_duplicate:
            liner.send('{}: {}\n'.format(
                readfile_field, _display_path(readfile)))

    # Return readfile validity
    return all([
        readfile_exists,
        readfile_is_fastq,
        not readfile_duplicate])

def get_optional_readfile_validity(
    readfile,
    readfile_field,
    paired_readfile,
    liner):
    '''
    Determine if an optional readfile exists,
    is non-empty and of FastQ type? If provided,
    check if given readfile is a duplicate of
    paired_readfile? Internal use only.

    :: readfile
       type - string
       desc - path to FastQ file storing
              reads
    :: readfile_field
       type - string
       desc - readfile fieldname used in
              printing
    :: paired_readfile
       type - string / None
       desc - path to paired FastQ file
              storing reads
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # readfile is None
    if readfile is None:
        liner.send(
            '{}: None Specified\n'.format(
                readfile_field))
        return True

    # Regular readfile Validation
    return get_readfile_validity(
        readfile=readfile,
        readfile_field=readfile_field,
        paired_readfile=paired_readfile,
        liner=liner)

def get_inzip_validity(
    inzip,
    inzip_suffix,
    inzip_field,
    inzip_type,
    liner):
    '''
    Determine if zipfile is valid.
    Internal use only.

    :: inzip
       type - string
       desc - path to compressed zipfile
    :: inzip_suffix
       type - string
       desc - zipfile suffix adjustment
    :: inzip_field
       type - string
       desc - zipfile fieldname used
              in printing
    :: inzip_type
       type - string
       desc - zipfile type name string
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    fname = _display_path(inzip)
    # Use space separator for numbered list items (field ends with ])
    sep = '' if inzip_field.rstrip().endswith(']') else ':'
    liner.send(
        '{}{} Loading {} ...'.format(
            inzip_field, sep, fname))

    # zipfile exists and non-empty?
    inzip_exists = get_infile_validity(
        infile=inzip,
        infile_suffix=inzip_suffix,
        infile_field=inzip_field,
        liner=liner)

    # inzip is zipfile?
    inzip_is_zipfile = False
    if inzip_exists:
        inzip = ut.get_adjusted_path(
            path=inzip,
            suffix=inzip_suffix)
        inzip_is_zipfile = zf.is_zipfile(
            filename=inzip)
        if not inzip_is_zipfile:
            liner.send(
                '{}{} {} [INVALID {} FILE]\n'.format(
                    inzip_field, sep, _display_path(inzip), inzip_type))

    # inzip is good?
    inzip_is_good = False
    archive = None
    if inzip_is_zipfile:
        archive = ut.get_archive(
            arcfile=inzip)
        inzip_is_good = True
        # inzip_is_good = archive.testzip() is None
        # if not inzip_is_good:
        #     liner.send(
        #         '{}: {} [CORRUPT {} FILE]\n'.format(
        #             inzip_field, inzip, inzip_type))

    # inzip valid?
    return (all([
        inzip_exists,
        inzip_is_zipfile,
        inzip_is_good]),
        inzip,
        archive)

def get_indexfile_validity(
    indexfile,
    indexfile_field,
    associated,
    liner):
    '''
    Determine if indexfile points to a
    valid zipfile with indexed objects.
    Internal use only.

    :: indexfile
       type - string
       desc - path to compressed zipfile
              storing prepared structures
              and data models
              (suffix='.oligopool.index')
    :: indexfile_field
       type - string
       desc - indexfile fieldname used in
              printing
    :: associated
       type - boolean
       desc - if True then indexfile is
              expected to have associates
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # indexfile exists and non-empty?
    (indexfile_ok_format,
    indexfile,
    archive) = get_inzip_validity(
        inzip=indexfile,
        inzip_suffix='.oligopool.index',
        inzip_field=indexfile_field,
        inzip_type='INDEX',
        liner=liner)

    # indexfile content ok?
    indexfile_ok_content = False
    # Use space separator for numbered list items (field ends with ])
    sep = '' if indexfile_field.rstrip().endswith(']') else ':'
    if indexfile_ok_format:
        try:
            indexed = set([
                'ID.map',
                'barcode.model',
                'meta.map'])
            if associated:
                indexed.add('associate.map')
            if not (indexed <= set(archive.namelist())):
                raise
            variantcount = ut.loaddict(
                archive=archive,
                dfile='meta.map')['variantcount']
        except:
            liner.send(
                '{}{} {} [INVALID INDEX FILE]\n'.format(
                    indexfile_field, sep, _display_path(indexfile)))
        else:
            liner.send(
                '{}{} {} w/ {:,} Variant(s)\n'.format(
                    indexfile_field, sep, _display_path(indexfile), variantcount))
            indexfile_ok_content = True
        finally:
            archive.close()

    # indexfile valid?
    return all([
        indexfile_ok_format,
        indexfile_ok_content])

def get_indexfiles_validity(
    indexfiles,
    indexfiles_field,
    associated,
    liner):
    '''
    Determine if all indexfiles are valid.
    Internal use only.

    :: indexfiles
       type - string / list
       desc - path to compressed zipfile
              storing prepared structures
              and data models
              (suffix='.oligopool.index')
    :: indexfiles_field
       type - string
       desc - indexfile fieldname used in
              printing
    :: associated
       type - boolean
       desc - if True then indexfiles are
              expected to have associates
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Do we have a single indexfile?
    if isinstance(indexfiles, str):
        return get_indexfile_validity(
            indexfile=indexfiles,
            indexfile_field=indexfiles_field,
            associated=associated,
            liner=liner)

    # Are indexfiles iterable?
    indexfiles_iterable = False
    indexfile_store = cx.deque()
    try:
        for indexfile in indexfiles:
            indexfile_store.append(indexfile)
        else:
            indexfiles_iterable = True
    except:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                indexfiles_field, indexfiles))

    # Are indexfiles valid?
    indexfiles_ok = True
    seen_indexfiles  = set()
    if indexfiles_iterable:

        # Prepare Spacing
        altspacing = ' '*len(indexfiles_field)

        # Show Header Update
        liner.send(
            '{}: {:,} Index Input(s)\n'.format(
                indexfiles_field, len(indexfile_store)))

        # Core Validation Loop
        idx = 0
        while indexfile_store:

            # Fetch indexfile
            indexfile = indexfile_store.popleft()
            idx += 1

            # Duplicate indexfile?
            if indexfile in seen_indexfiles:
                # Show Update
                liner.send(
                    '{}:  [{}] {} [DUPLICATE INDEX FILE]\n'.format(
                        altspacing, idx, _display_path(ut.get_adjusted_path(
                            path=indexfile,
                            suffix='.oligopool.index'))))
                # Update Global Validity
                indexfiles_ok = indexfiles_ok and False
                # Next!
                continue

            # Record indexfile
            else:
                seen_indexfiles.add(indexfile)

            # Get Current Validity
            idxfile_valid = get_indexfile_validity(
                indexfile=indexfile,
                indexfile_field='{}:  [{}]'.format(altspacing, idx),
                associated=associated,
                liner=liner)

            # Update Global Validity
            indexfiles_ok = indexfiles_ok and idxfile_valid

    # Return Results
    return indexfiles_iterable and indexfiles_ok

def get_parsed_packfile(
    packfile,
    packfile_field,
    liner):
    '''
    Determine if packfile points to a
    valid zipfile with .pack files.
    Also, return packcount.
    Internal use only.

    :: packfile
       type - string
       desc - path to compressed zipfile
              storing read packs
              (suffix='.oligopool.pack')
    :: packfile_field
       type - string
       desc - packfile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # packfile exists and non-empty?
    (packfile_ok_format,
    packfile,
    archive) = get_inzip_validity(
        inzip=packfile,
        inzip_suffix='.oligopool.pack',
        inzip_field=packfile_field,
        inzip_type='PACK',
        liner=liner)

    # packfile content ok and non-empty?
    packfile_ok_content = False
    packcount = 0
    packfile_nonempty = False
    if packfile_ok_format:
        try:

            for cpath in archive.namelist():
                if not (cpath.endswith('.pack') or \
                        cpath.endswith('.stat')):
                    raise
                packcount += 1

            packcount -= 1 # Stat File isn't Pack

            packstat = ut.loaddict(
                archive=archive,
                dfile='packing.stat')

            if packcount != packstat['pack_count']:
                raise

            packfile_nonempty = packcount > 0

        except:
            # Invalid / corrupt packfile
            liner.send(
                '{}: {} [INVALID PACK FILE]\n'.format(
                    packfile_field, _display_path(packfile)))
        else:
            # Empty packfile
            if not packfile_nonempty:
                liner.send(
                    '{}: {} w/ {} Read Packs [EMPTY PACK FILE]\n'.format(
                        packfile_field, _display_path(packfile), packcount))
            # Packfile has packs
            else:
                liner.send(
                    '{}: {} w/ {} Read Packs\n'.format(
                        packfile_field, _display_path(packfile), packcount))
                packfile_ok_content = True
        finally:
            archive.close()

    # packfile valid?
    return (all([
        packfile_ok_format,
        packfile_ok_content,
        packfile_nonempty]),
        packcount)

def get_parsed_data_info(
    data,
    data_field,
    required_fields,
    liner,
    allow_missing_cols=None,
    id_aliases=None):
    '''
    Determine if given data is a valid,
    non-empty CSV file or a DataFrame.
    Internal use only.

    :: data
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing information
    :: data_field
       type - string
       desc - data fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: liner
       type - coroutine
       desc - dynamic printing
    :: allow_missing_cols
       type - set / None
       desc - column names allowed to
              contain missing values
    '''

    data_field = _normalize_field(data_field)

    # data flags
    data_type   = None
    data_is_df  = False
    data_is_csv = False

    # df flags
    df = None
    df_nonempty = False

    # data is a DataFrame?
    if isinstance(data, pd.DataFrame):
        df = data.copy().reset_index(
            drop=data.index.name is None)
        data_is_df = True

    # data is a valid CSV file?
    else:

        # data file exists and non-empty?
        data_is_file = get_infile_validity(
            infile=data,
            infile_suffix=None,
            infile_field=data_field,
            liner=liner)

        # data is a CSV file?
        if data_is_file:

            # CSV file read attempt
            try:
                # Parse full df
                df = pd.read_csv(
                    filepath_or_buffer=data,
                    sep=',',
                    header=0,
                    engine='c')

            # Read and/or Indexing Unsuccesful
            except:
                liner.send('{}: {} [INVALID CSV FILE]\n'.format(
                    data_field,
                    _display_path(data)))

            # Read succesful
            else:
                data_is_csv = True

    # df was extracted?
    df_extracted = data_is_df or data_is_csv

    # Allow alternate ID column names (e.g., DegenerateID) by renaming to 'ID'
    # before running the standard missing-value and indexing checks.
    if df_extracted and \
       not df is None and \
       'ID' not in df.columns and \
       id_aliases:
        for alias in id_aliases:
            if alias in df.columns:
                df.rename(columns={alias: 'ID'}, inplace=True)
                break

    # Update data and data_type
    if   data_is_df:
        data      = 'DataFrame'
        data_type = 'DATAFRAME'
    elif data_is_csv:
        data      = _display_path(data)
        data_type = 'CSV FILE'

    # df is non-empty?
    df_nonempty = False

    if df_extracted:

        # Compute emptiness
        df_nonempty = not df.empty and \
                      len(df.columns) > 1

        if not df_nonempty:
            liner.send(
                '{}: {} w/ {:,} Record(s) [{} IS EMPTY]\n'.format(
                    data_field,
                    data,
                    len(df.index),
                    data_type))

    # df columns have missing values?
    df_no_missing_vals = False

    if df_nonempty:

        # Check column-wise missing values
        for col in df.columns:
            if allow_missing_cols and col in allow_missing_cols:
                continue
            missing_mask = ut.get_missing_mask(
                series=df[col],
                allow_dash=False)
            if missing_mask.any():
                examples = ut.get_row_examples(
                    df=df,
                    invalid_mask=missing_mask,
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '{}: {} w/ {:,} Record(s) [MISSING VALUES IN COLUMN=\'{}\']{}\n'.format(
                        data_field,
                        data,
                        len(df.index),
                        col,
                        example_note))
                break

        # No missing values
        else:
            df_no_missing_vals = True

    # df indexible?
    df_indexible = False

    if df_no_missing_vals:

        # Try indexing on unique ID
        try:

            # Reset and (Re-)Index by ID
            # print(df.index.is_unique)
            df.set_index(
                keys='ID',
                inplace=True)

        # Indexing unsuccessful
        except Exception:

            # Unindexible df
            liner.send(
                '{}: {} w/ {:,} Record(s) [INVALID OR MISSING COLUMN=\'ID\']\n'.format(
                    data_field,
                    data,
                    len(df.index)))

            # Indexing failed
            df_indexible = False

        else:

            # Assert ID keys are unique
            if not df.index.is_unique:
                duplicates = []
                for value in df.index[df.index.duplicated()]:
                    sval = str(value)
                    if sval not in duplicates:
                        duplicates.append(sval)
                    if len(duplicates) >= 5:
                        break
                example_note = ut.format_row_examples(
                    duplicates,
                    label='Duplicate ID examples')
                liner.send(
                    '{}: {} w/ {:,} Record(s) [NON-UNIQUE COLUMN=\'ID\']{}\n'.format(
                        data_field,
                        data,
                        len(df.index),
                        example_note))
                df_indexible = False
            else:
                # Normalize ID values to strings for consistent downstream behavior.
                df.index = df.index.map(str)

                # Ensure uniqueness is preserved after casting (e.g., 1 and "1").
                if not df.index.is_unique:
                    duplicates = []
                    for value in df.index[df.index.duplicated()]:
                        sval = str(value)
                        if sval not in duplicates:
                            duplicates.append(sval)
                        if len(duplicates) >= 5:
                            break
                    example_note = ut.format_row_examples(
                        duplicates,
                        label='Duplicate ID examples')
                    liner.send(
                        '{}: {} w/ {:,} Record(s) [NON-UNIQUE COLUMN=\'ID\' AFTER STRING CAST]{}\n'.format(
                            data_field,
                            data,
                            len(df.index),
                            example_note))
                    df_indexible = False
                else:
                    # Everything checked out
                    df_indexible = True

    # df contains required columns?
    df_contains_required_cols = False

    if df_indexible:

        # Required fields present?
        if not required_fields is None:

            # Track requirement fulfillment
            required_fields_found = 0

            # Get df columns
            present_cols = set([c.lower() for c in df.columns])
            present_cols.update((df.index.name.lower(),))

            # Loop through required fields
            for required_field in required_fields:

                # Required field in present df columns
                if required_field.lower() in present_cols:
                    required_fields_found += 1
                # Missing field!
                else:
                    break

            # Requirements fulfilled
            if required_fields_found == len(required_fields):
                df_contains_required_cols = True

            # Requirement failed
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [MISSING COLUMN=\'{}\']\n'.format(
                        data_field,
                        data,
                        len(df.index),
                        required_field))

        # No requirements specified
        else:
            df_contains_required_cols = True

    # Compute final validity
    df_valid = all([
        df_extracted,
        df_nonempty,
        df_no_missing_vals,
        df_indexible,
        df_contains_required_cols])

    # Return Results
    data_name = data
    return (df, data_name, df_valid)

def get_parsed_verify_indata_info(
    data,
    data_field,
    liner):
    '''
    Parse and validate input data for verify module.
    Checks for ID column and at least one DNA column (ATGC only).
    Internal use only.

    :: data
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing information
    :: data_field
       type - string
       desc - data fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing

    Returns:
        (df, data_name, dna_columns, valid)
        - df: parsed DataFrame (ID as index) or None
        - data_name: 'DataFrame' or filename
        - dna_columns: list of DNA column names found
        - valid: True if validation passed
    '''

    # First pass: basic validation via get_parsed_data_info
    (df,
    data_name,
    base_valid) = get_parsed_data_info(
        data=data,
        data_field=data_field,
        required_fields=('ID',),
        liner=liner)

    # If base validation failed, return early
    if not base_valid:
        return (df, data_name, [], False)

    # Detect DNA columns
    dna_columns = []
    used_complete_oligo = False

    # Check for CompleteOligo first
    if 'CompleteOligo' in df.columns:
        # Validate CompleteOligo contains DNA
        series = df['CompleteOligo']
        if series.map(lambda x: isinstance(x, str)).all():
            upper = series.str.upper()
            is_dna = upper.map(lambda x: ut.is_DNA(seq=x, dna_alpha=ut.dna_alpha))
            if is_dna.all():
                dna_columns = ['CompleteOligo']
                used_complete_oligo = True

    # If no CompleteOligo, scan all columns
    if not used_complete_oligo:
        for col in df.columns:
            series = df[col]
            # Skip non-string columns
            if not series.map(lambda x: isinstance(x, str)).all():
                continue
            # Check if all values are pure ATGC (plus gaps)
            upper = series.str.upper()
            is_dna = upper.map(lambda x: ut.is_DNA(seq=x, dna_alpha=ut.dna_alpha))
            if is_dna.all():
                dna_columns.append(col)

    # Validate at least one DNA column found
    if not dna_columns:
        liner.send(
            '{}: {} w/ {:,} Record(s) [NO DNA COLUMNS DETECTED]\n'.format(
                data_field,
                data_name,
                len(df.index)))
        return (df, data_name, [], False)

    # All validation passed
    return (df, data_name, dna_columns, True)

def get_parsed_indata_info(
    indata,
    indata_field,
    required_fields,
    precheck,
    liner,
    allow_missing_cols=None):
    '''
    Determine if indata consisting of DNA columns
    only is valid. Internal use only.

    :: indata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing information
    :: indata_field
       type - string
       desc - indata fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: precheck
       type - boolean
       desc - if False prints content description
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    :: allow_missing_cols
       type - set / None
       desc - column names allowed to
              contain missing values
    '''

    indata_field = _normalize_field(indata_field)

    # Is indata valid CSV or DataFrame?
    (df,
    data_name,
    df_valid) = get_parsed_data_info(
        data=indata,
        data_field=indata_field,
        required_fields=required_fields,
        liner=liner,
        allow_missing_cols=allow_missing_cols)

    # df columns contain DNA strings only?
    df_contains_DNA_only = False

    if df_valid:

        # Are all entries DNA strings?

        non_DNA_found = False

        # Loop through all columns
        for column in df.columns:

            # Missing allowed?
            missing_mask = None
            if allow_missing_cols and column in allow_missing_cols:
                missing_mask = ut.get_missing_mask(
                    series=df[column],
                    allow_dash=True)

            # Loop through all entires in column
            for idx, value in enumerate(df[column]):

                # Skip missing in allowed column
                if missing_mask is not None and missing_mask[idx]:
                    continue

                # A Non-DNA entry?
                if not ut.is_DNA(seq=value):

                    # Set flag
                    non_DNA_found = True

                    # Terminate inner loop
                    break

            # We encountered Non-DNA entry
            if non_DNA_found:
                # Terminate outer loop
                break

        df_contains_DNA_only = not non_DNA_found

        # SHow update
        if not df_contains_DNA_only:
            invalid_mask = ~df[column].map(
                lambda x: ut.is_DNA(seq=x))
            examples = ut.get_row_examples(
                df=df,
                invalid_mask=invalid_mask,
                id_col='ID',
                limit=5)
            example_note = ut.format_row_examples(examples)
            liner.send(
                '{}: {} w/ {:,} Record(s) [NON-DNA VALUE=\'{}\' IN COLUMN=\'{}\']{}\n'.format(
                    indata_field,
                    data_name,
                    len(df.index),
                    value,
                    column,
                    example_note))

    # Compute final validity
    df_valid = all([
        df_valid,
        df_contains_DNA_only])

    # Is df valid?
    if df_valid:

        # Uppercase all DNA Strings
        for col in df.columns:
            df[col] = df[col].str.upper()

        # Show update?
        if not precheck:
            liner.send(
                '{}: {} w/ {:,} Record(s)\n'.format(
                    indata_field,
                    data_name,
                    len(df.index)))

    else:
        # Erase df
        df = None

    # Return data validity
    if precheck:
        return (df, data_name, df_valid)
    return (df, df_valid)

def get_outfile_validity(
    outfile,
    outfile_suffix,
    outfile_field,
    liner):
    '''
    Determine if outfile points to an existing
    empty file or non-existent but creatable path?
    Internal use only.

    :: outfile
       type - string
       desc - an output file storing
              computed information
    :: outfile_suffix
       type - string
       desc - required outfile suffix
    :: outfile_field
       type - string
       desc - outfile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    outfile_field = _normalize_field(outfile_field)

    # Append suffix to path?

    # outfile is an existing empty file or
    # non-existent but creatable path?
    outfile_status = ut.get_path_status(
        path=outfile,
        suffix=outfile_suffix,
        readable=False,
        writable=True,
        creatable=True)

    outfile_valid = False

    # outfile is non-string type?
    if outfile_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                outfile_field, _display_path(outfile)))

    else:

        # Adjust outfile with suffix
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix=outfile_suffix)

        # outfile is invalid
        if outfile_status in (2, 11):
            liner.send(
                '{}: {} [WRITE PERMISSION DENIED]\n'.format(
                    outfile_field, _display_path(outfile)))

        elif outfile_status == 4:
            liner.send(
                '{}: {} [FILE ALREADY EXISTS]\n'.format(
                    outfile_field, _display_path(outfile)))

        elif outfile_status == 'X':
            liner.send(
                '{}: {} [FILE IS SPECIAL]\n'.format(
                    outfile_field, _display_path(outfile)))

        elif 6 <= outfile_status <= 8:
            liner.send(
                '{}: {} [FILE IS DIRECTORY]\n'.format(
                    outfile_field, _display_path(outfile)))

        # outfile is valid
        elif outfile_status in (3, 10):
            liner.send('{}: {}\n'.format(
                outfile_field, _display_path(outfile)))
            outfile_valid = True

    # Return outfile validity
    return outfile_valid

def get_optional_outfile_validity(
    outfile,
    outfile_suffix,
    outfile_field,
    liner):
    '''
    Determine if an optional outfile points to an existing
    empty file or non-existent but creatable path?
    Internal use only.

    :: outfile
       type - string / None
       desc - an optional output file storing
              computed information
    :: outfile_suffix
       type - string
       desc - required outfile suffix
    :: outfile_field
       type - string
       desc - outfile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    outfile_field = _normalize_field(outfile_field)

    # outfile is None - feature disabled
    if outfile is None:
        liner.send('{}: Disabled\n'.format(outfile_field))
        return True

    # Regular outfile validation
    return get_outfile_validity(
        outfile=outfile,
        outfile_suffix=outfile_suffix,
        outfile_field=outfile_field,
        liner=liner)

def get_outdir_validity(
    outdir,
    outdir_suffix,
    outdir_field,
    liner):
    '''
    Determine if outdir points to an existing
    empty directory or non-existent but creatable
    path? Internal use only.

    :: outdir
       type - string
       desc - an output directory storing
              computed information
    :: outdir_suffix
       type - string
       desc - required outdir suffix
    :: outdir_field
       type - string
       desc - outdir fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    outdir_field = _normalize_field(outdir_field)

    # Append suffix to path?

    # outdir is an existing empty directory
    # or non-existent but creatable path?
    outdir_status = ut.get_path_status(
        path=outdir,
        suffix=outdir_suffix,
        readable=False,
        writable=True,
        creatable=True)

    outdir_valid = False

    # outdir is non-string type?
    if outdir_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                outdir_field, _display_path(outdir)))

    else:

        # Adjust outdir with suffix
        outdir = ut.get_adjusted_path(
            path=outdir,
            suffix=outdir_suffix)

        # outdir is invalid
        if outdir_status in (6, 11):
            liner.send(
                '{}: {} [WRITE PERMISSION DENIED]\n'.format(
                    outdir_field, _display_path(outdir)))

        elif outdir_status == 8:
            liner.send(
                '{}: {} [DIRECTORY ALREADY EXISTS]\n'.format(
                    outdir_field, _display_path(outdir)))

        elif outdir_status == 'X':
            liner.send(
                '{}: {} [DIRECTORY IS SPECIAL]\n'.format(
                    outdir_field, _display_path(outdir)))

        elif 1 <= outdir_status <= 4:
            liner.send(
                '{}: {} [DIRECTORY IS FILE]\n'.format(
                    outdir_field, _display_path(outdir)))

        # outdir is valid
        elif outdir_status in (7, 10):
            liner.send('{}: {}\n'.format(
                outdir_field, _display_path(outdir)))
            outdir_valid = True

    # Return outdir validity
    return outdir_valid

def get_outdf_validity(
    outdf,
    outdf_suffix,
    outdf_field,
    liner):
    '''
    Determine if outdf points to a valid outfile
    if specified. Internal use only.

    :: outdf
       type - string / None
       desc - output file storing
              computed information
    :: outdf_suffix
       type - string
       desc - required outdf suffix
    :: outfile_field
       type - string
       desc - dipath fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    outdf_field = _normalize_field(outdf_field)

    # outdf is None?
    if outdf is None:
        liner.send(
            '{}: In-Memory DataFrame\n'.format(
                outdf_field))
        return True

    # outdf is valid outfile?
    return get_outfile_validity(
        outfile=outdf,
        outfile_suffix=outdf_suffix,
        outfile_field=outdf_field,
        liner=liner)

def get_parsed_flag_info(
    flag,
    flag_field,
    flag_desc_off,
    flag_desc_on,
    liner,
    precheck=False,
    flag_valid=None,
    flag_on=None):
    '''
    Determine if a given input is a valid boolean-like flag (0/1/True/False) and optionally
    print a standardized line to the module argument recap. Internal use only.

    :: flag
       type - bool / int / np.bool_
       desc - the flag input
    :: flag_field
       type - string
       desc - fieldname used in printing
    :: flag_desc_off
       type - string
       desc - description printed when flag is False
    :: flag_desc_on
       type - string
       desc - description printed when flag is True
    :: liner
       type - coroutine
       desc - dynamic printing
    :: precheck
       type - bool
       desc - if True, skip printing (default: False)
    :: flag_valid
       type - bool / None
       desc - optional precomputed validity
    :: flag_on
       type - bool / None
       desc - optional precomputed boolean value
    '''

    flag_field = _normalize_field(flag_field)

    # Compute validity if not supplied
    if (flag_valid is None) or (flag_on is None):
        flag_valid = isinstance(flag, (bool, int, np.bool_)) and \
            (flag in (0, 1, True, False))
        flag_on = bool(flag) if flag_valid else False

    # Print if requested
    if not precheck:
        if flag_valid:
            liner.send('{}: {}\n'.format(
                flag_field,
                [flag_desc_off, flag_desc_on][flag_on]))
        else:
            liner.send('{}: {} [INPUT TYPE IS INVALID]\n'.format(
                flag_field,
                flag))

    return (flag_on, flag_valid)

def get_parsed_random_seed_info(
    random_seed,
    random_seed_field,
    liner):
    '''
    Determine if random_seed is a valid RNG seed (non-negative int / None) and print a standardized
    line to the module argument recap. Internal use only.

    :: random_seed
       type - int / np.integer / None
       desc - RNG seed (None means unseeded RNG)
    :: random_seed_field
       type - string
       desc - fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    random_seed_field = _normalize_field(random_seed_field)

    # No seed specified
    if random_seed is None:
        liner.send('{}: None Specified\n'.format(random_seed_field))
        return (None, True)

    # Reject booleans explicitly (bool is a subclass of int)
    if isinstance(random_seed, (bool, np.bool_)):
        liner.send('{}: {} [MUST BE INTEGER OR NONE]\n'.format(
            random_seed_field, random_seed))
        return (random_seed, False)

    # Accept Python / numpy integers
    if isinstance(random_seed, (int, np.integer)):
        random_seed = int(random_seed)
        if random_seed < 0:
            liner.send('{}: {} [MUST BE NON-NEGATIVE]\n'.format(
                random_seed_field, random_seed))
            return (random_seed, False)
        liner.send('{}: {:,}\n'.format(
            random_seed_field, random_seed))
        return (random_seed, True)

    # Invalid type
    liner.send('{}: {} [MUST BE INTEGER OR NONE]\n'.format(
        random_seed_field, random_seed))
    return (random_seed, False)

def get_parsed_column_info(
    col,
    df,
    col_field,
    col_desc,
    col_type,
    adjcol,
    adjval,
    iscontext,
    typecontext,
    liner,
    allow_existing=False):
    '''
    Determine if col is a valid column in/for df
    depending on col_type. Internal use only.

    :: col
       type - string / None
       desc - column name to parse
    :: df
       type - pd.DataFrame
       desc - DataFrame to store output in
              or get input from
    :: col_field
       type - string
       desc - col fieldname used in
              printing
    :: col_desc
       type - string
       desc - col description used in
              printing
    :: col_type
       type - string
       desc - if 0, treat the column as input,
              otherwise, treat as output column
    :: adjcol
       type - string / None
       desc - ensure col is adjacent to adjcol
    :: adjval
       type - integer / None
       desc - ensure adjcol adjacent by adjval
    :: iscontext
       type - boolean
       desc - if True treat column as context,
              otherwise, column is a singular
              and context free element
    :: typecontext
       type - integer / None
       desc - if 0, treat context as left,
              otherwse, treat as right
    :: liner
       type - coroutine
       desc - dynamic printing
    :: allow_existing
       type - boolean
       desc - if True, allow output
              column to already exist
    '''

    col_field = _normalize_field(col_field)

    # Is col None?
    if col is None:
        # OK for Input Column
        if col_type == 0:
            liner.send(
                '{}: None Specified\n'.format(
                    col_field,
                    col_desc,
                    col))
            return (None, True)

    # Is col string?
    col_is_string = False

    if isinstance(col, str):
        col_is_string = True

    else:
        liner.send(
            '{}: {} \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                col_field,
                col_desc,
                col))

    # Is col non-existent?
    col_existence_valid = False

    if col_is_string:

        # Does df exist? No
        if not isinstance(df, pd.DataFrame):

            # There's no df to check col in
            liner.send(
                '{}: {} \'{}\' [COLUMN EXISTENCE UNKNOWN]\n'.format(
                    col_field,
                    col_desc,
                    col))

            # Guess, that means we don't know
            # the status of column specified
            col_existence_valid = False

        # Does df exist? Yes
        else:

            # Get col existence
            (col_existence,
            col_idx) = ut.get_col_exist_idx(
                col=col,
                df=df)

            # Is col_existence valid?
            if  ((col_type == 0 and    #  Input Column
                  col_existence) or    #  Input Column
                 (col_type == 1 and    # Output Column
                  (not col_existence or allow_existing))): # Output Column
                col_existence_valid = True
            else:
                col_existence_valid = False

            # Show update for invalidation
            if not col_existence_valid:

                # Input Column
                if col_type == 0:
                    liner.send(
                    '{}: {} \'{}\' [COLUMN DOES NOT EXIST]\n'.format(
                        col_field,
                        col_desc,
                        col))

                # Output Column
                else:
                    liner.send(
                    '{}: {} \'{}\' [COLUMN ALREADY EXISTS]\n'.format(
                        col_field,
                        col_desc,
                        col))

    # Is col adjacent?
    col_is_adjacent = False

    if col_type == 0       and \
       col_existence_valid and \
       not adjcol is None  and \
       not adjval is None:

        # Get adjcol existence and index
        (adj_existence,
        adj_idx) = ut.get_col_exist_idx(
            col=adjcol,
            df=df)

        # adjcol non-existent?
        if not adj_existence:
            liner.send(
                '{}: {} \'{}\' [ADJACENT COLUMN \'{}\' DOES NOT EXIST]\n'.format(
                    col_field,
                    col_desc,
                    col,
                    adjcol))

        # adjcol more than adjval away?
        elif adj_idx - col_idx != adjval:
            liner.send(
                '{}: {} \'{}\' [COLUMN \'{}\' NON-ADJACENT]\n'.format(
                    col_field,
                    col_desc,
                    col,
                    adjcol))

        # All conditions met!
        else:
            col_is_adjacent = True

    else:
        col_is_adjacent = True


    # Is col valid?
    col_valid = all([col_is_string,
        col_existence_valid,
        col_is_adjacent])

    if col_valid:
        if col_type == 1 and allow_existing:
            # Output Column (may already exist)
            suffix = ' (Existing)' if col in df.columns else ''
            liner.send(
                '{}: {} \'{}\'{}\n'.format(
                    col_field,
                    col_desc,
                    col,
                    suffix))
        else:
            liner.send(
                '{}: {} \'{}\'\n'.format(
                    col_field,
                    col_desc,
                    col))

        # Output Column
        if col_type == 1:
            parsedcol = None

        # Input Column
        else:
            # Non-Context Extraction
            if not iscontext:
                parsedcol = df[col].str.upper()
            # Context Extraction
            else:
                if typecontext == 0: #  Left Context
                    parsedcol = df.loc[:, :col]
                else:                # Right Context
                    parsedcol = df.loc[:, col:]
                # Robust to allowed missing values in output columns (Patch Mode fills):
                # treat NaN/None as empty during context concatenation.
                parsedcol = parsedcol.fillna('')
                parsedcol = parsedcol.sum(
                    axis=1).str.upper().str.replace(
                        '-', '')
    else:
        parsedcol = None

    # Return based on col_type
    if col_type == 1:
        return col_valid
    else:
        return (parsedcol,
            col_valid)

def get_parsed_cross_barcode_info(
    crosscols,
    crosscols_field,
    mindist,
    mindist_field,
    barcodelen,
    outcol,
    df,
    liner):
    '''
    Determine if cross barcode constraints are valid.
    Internal use only.

    :: crosscols
       type - list / None
       desc - list of column names in df storing existing barcodes
    :: crosscols_field
       type - string
       desc - crosscols fieldname used in printing
    :: mindist
       type - integer / None
       desc - minimum cross-set Hamming distance
    :: mindist_field
       type - string
       desc - mindist fieldname used in printing
    :: barcodelen
       type - integer
       desc - barcode length
    :: outcol
       type - string
       desc - output barcode column name
    :: df
       type - pd.DataFrame / None
       desc - input pandas DataFrame
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    crosscols_field = _normalize_field(crosscols_field)
    mindist_field = _normalize_field(mindist_field)

    # For barcodes, accept only strict ATGC (no gaps / no degeneracy).
    barcode_alpha = set('ATGC')

    # Normalize crosscols input
    if isinstance(crosscols, str):
        crosscols = [crosscols]

    # Cross mode disabled
    if crosscols is None and mindist is None:
        liner.send(
            '{}: None Specified\n'.format(
                crosscols_field))
        liner.send(
            '{}: None Specified\n'.format(
                mindist_field))
        return (None, None, True)

    # Symmetry check
    if (crosscols is None) ^ (mindist is None):
        if crosscols is None:
            liner.send(
                '{}: None Specified [MISSING CROSS BARCODE COLUMNS]\n'.format(
                    crosscols_field))
        else:
            if not isinstance(crosscols, (list, tuple)):
                liner.send(
                    '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                        crosscols_field,
                        crosscols))
            else:
                liner.send(
                    '{}: Input from Column(s) \'{}\'\n'.format(
                        crosscols_field,
                        ', '.join(str(c) for c in crosscols)))
        if mindist is None:
            liner.send(
                '{}: None Specified [MISSING MINIMUM CROSS DISTANCE]\n'.format(
                    mindist_field))
        else:
            if isinstance(mindist, nu.Integral):
                liner.send(
                    '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
                        mindist_field,
                        mindist))
            else:
                liner.send(
                    '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                        mindist_field,
                        mindist))
        return (None, None, False)

    # Type check: crosscols must be a list/tuple of strings
    if not isinstance(crosscols, (list, tuple)) or not crosscols:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                crosscols_field,
                crosscols))
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                mindist_field,
                mindist))
        return (None, None, False)

    if not all(isinstance(col, str) for col in crosscols):
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                crosscols_field,
                crosscols))
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                mindist_field,
                mindist))
        return (None, None, False)

    crosscols = list(crosscols)

    # Output column conflict?
    if outcol in crosscols:
        liner.send(
            '{}: Input from Column(s) \'{}\' [COLUMN NAME CONFLICT=\'{}\']\n'.format(
                crosscols_field,
                ', '.join(crosscols),
                outcol))
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                mindist_field,
                mindist))
        return (None, None, False)

    # mindist type check
    if not isinstance(mindist, nu.Integral):
        liner.send(
            '{}: Input from Column(s) \'{}\'\n'.format(
                crosscols_field,
                ', '.join(crosscols)))
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                mindist_field,
                mindist))
        return (None, None, False)

    mindist = int(mindist)

    # mindist range check
    if mindist < 1 or mindist > barcodelen:
        liner.send(
            '{}: Input from Column(s) \'{}\'\n'.format(
                crosscols_field,
                ', '.join(crosscols)))
        liner.send(
            '{}: {} [VALUE OUT OF RANGE]\n'.format(
                mindist_field,
                mindist))
        return (None, None, False)

    # df exists?
    if not isinstance(df, pd.DataFrame):
        liner.send(
            '{}: Input from Column(s) \'{}\' [COLUMN EXISTENCE UNKNOWN]\n'.format(
                crosscols_field,
                ', '.join(crosscols)))
        liner.send(
            '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
                mindist_field,
                mindist))
        return (crosscols, mindist, False)

    # Column existence and content validity
    for col in crosscols:

        if col not in df.columns:
            liner.send(
                '{}: Input from Column(s) \'{}\' [MISSING COLUMN=\'{}\']\n'.format(
                    crosscols_field,
                    ', '.join(crosscols),
                    col))
            liner.send(
                '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
                    mindist_field,
                    mindist))
            return (None, None, False)

        uniques = ut.get_uniques(
            iterable=df[col],
            typer=tuple)

        for value in uniques:
            if not isinstance(value, str) or not ut.is_DNA(
                seq=value,
                dna_alpha=barcode_alpha):
                liner.send(
                    '{}: Input from Column(s) \'{}\' [NON-DNA VALUE IN COLUMN=\'{}\']\n'.format(
                        crosscols_field,
                        ', '.join(crosscols),
                        col))
                liner.send(
                    '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
                        mindist_field,
                        mindist))
                return (None, None, False)
            if len(value) != barcodelen:
                liner.send(
                    '{}: Input from Column(s) \'{}\' [INVALID LENGTH IN COLUMN=\'{}\']\n'.format(
                        crosscols_field,
                        ', '.join(crosscols),
                        col))
                liner.send(
                    '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
                        mindist_field,
                        mindist))
                return (None, None, False)

    # Show updates
    liner.send(
        '{}: Input from Column(s) \'{}\'\n'.format(
            crosscols_field,
            ', '.join(crosscols)))
    liner.send(
        '{}: At least {} Mismatch(es) per Cross-Set Barcode Pair\n'.format(
            mindist_field,
            mindist))

    return (crosscols, mindist, True)

def get_parsed_exseqs_info(
    exseqs,
    exseqs_field,
    exseqs_desc,
    df_field,
    required,
    liner):
    '''
    Determine if given excluded sequences are valid.
    Internal use only.

    :: exseqs
       type - iterable / string / pd.DataFrame / None
       desc - iterable of DNA strings to be excluded
              from barcodes and around edges; optionally
              a DataFrame or a path to a CSV file with
              such sequences
    :: exseqs_field
       type - string
       desc - exseqs fieldname used in printing
    :: exseqs_desc
       type - string
       desc - exseqs description used in printing
    :: df_field
       type - string
       desc - name of the column in DataFrame storing
              excluded sequences
    :: required
       type - boolean
       desc - if True then exseqs cannot be None
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    exseqs_field = _normalize_field(exseqs_field)

    # Is exseqs None?
    if exseqs is None:
        liner.send(
            '{}: 0 {}\n'.format(
                exseqs_field,
                exseqs_desc,))
        # Required but missing
        if required:
            return (None, False)
        # Non required so OK
        else:
            return (None, True)

    # Is exseqs string?
    if isinstance(exseqs, str):

        # Is exseqs a single DNA string?
        if ut.is_DNA(seq=exseqs.upper()):
            liner.send(
                '{}: 1 {}\n'.format(
                    exseqs_field,
                    exseqs_desc,))
            return ([exseqs.upper()], True)

        # Is exseqs a FASTA file?
        fasta_exts = ('.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz')
        is_fasta = exseqs.lower().endswith(fasta_exts)

        # Also check content if extension doesn't match
        if not is_fasta:
            import os
            if os.path.isfile(exseqs):
                try:
                    import gzip
                    opener = gzip.open if exseqs.endswith('.gz') else open
                    with opener(exseqs, 'rt') as f:
                        first_char = f.read(1)
                        is_fasta = (first_char == '>')
                except:
                    pass

        if is_fasta:
            try:
                import os
                import pyfastx
                if not os.path.isfile(exseqs):
                    liner.send(
                        '{}: {} [FILE NOT FOUND]\n'.format(
                            exseqs_field,
                            _display_path(exseqs)))
                    return (None, False)

                # Parse FASTA file using pyfastx (handles gzip automatically)
                sequences = [seq.upper() for _, seq in pyfastx.Fasta(exseqs, build_index=False)]

                # Get unique sequences
                exseqs = list(ut.get_uniques(iterable=sequences, typer=list))

                # Validate all are DNA strings
                for seq in exseqs:
                    if not ut.is_DNA(seq=seq):
                        liner.send(
                            '{}: {:,} {} [NON-DNA SEQUENCE IN FASTA]\n'.format(
                                exseqs_field,
                                len(exseqs),
                                exseqs_desc))
                        return (None, False)

                liner.send(
                    '{}: {:,} {} (from FASTA)\n'.format(
                        exseqs_field,
                        len(exseqs),
                        exseqs_desc))
                return (exseqs, True)

            except Exception as e:
                liner.send(
                    '{}: {} [INVALID FASTA FILE]\n'.format(
                        exseqs_field,
                        _display_path(exseqs)))
                return (None, False)

    # Is exseqs iterable?
    if not isinstance(exseqs, pd.DataFrame) and \
       not isinstance(exseqs, str):

        # Try extracting exseqs
        try:
            exseqs = list(map(lambda x: x.upper(), ut.get_uniques(
                iterable=exseqs,
                typer=list)))

        # Error during extraction
        except:
            liner.send(
                '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                    exseqs_field,
                    exseqs))
            return (None, False)

        # Ensure all motifs are DNA strings
        for seq in exseqs:

            # Non-DNA element found!
            if not ut.is_DNA(seq=seq):
                liner.send(
                    '{}: {:,} {} [NON-DNA SEQUENCE=\'{}\']\n'.format(
                        exseqs_field,
                        len(exseqs),
                        exseqs_desc,
                        seq))
                return (None, False)

        # No error in extraction
        # All sequences are DNA strings
        liner.send(
            '{}: {:,} {}\n'.format(
                exseqs_field,
                len(exseqs),
                exseqs_desc))

        return (exseqs, True)

    # Is exseqs a CSV file or DataFrame?
    df = None
    df_valid = False

    # Handle DataFrame directly
    if isinstance(exseqs, pd.DataFrame):
        df = exseqs.copy()
        data_name = 'DataFrame'
    # Handle CSV file path
    elif isinstance(exseqs, str):
        try:
            df = pd.read_csv(exseqs, sep=',', header=0, engine='c')
            data_name = _display_path(exseqs)
        except:
            liner.send(
                '{}: {} [INVALID CSV FILE]\n'.format(
                    exseqs_field,
                    _display_path(exseqs)))
            return (None, False)

    # Check if df_field column exists
    if df is not None:
        if df_field not in df.columns:
            liner.send(
                '{}: {} [MISSING COLUMN=\'{}\']\n'.format(
                    exseqs_field,
                    data_name,
                    df_field))
            return (None, False)

        # Check for missing values in the motif column
        missing_mask = ut.get_missing_mask(series=df[df_field], allow_dash=False)
        if missing_mask.any():
            liner.send(
                '{}: {} w/ {:,} Record(s) [MISSING VALUES IN COLUMN=\'{}\']\n'.format(
                    exseqs_field,
                    data_name,
                    len(df.index),
                    df_field))
            return (None, False)

        # Extract unique motifs
        exseqs = list(map(lambda x: x.upper(), ut.get_uniques(
            iterable=df[df_field].to_list(),
            typer=list)))

        # Validate all are DNA strings
        for seq in exseqs:
            if not ut.is_DNA(seq=seq):
                liner.send(
                    '{}: {:,} {} [NON-DNA SEQUENCE=\'{}\']\n'.format(
                        exseqs_field,
                        len(exseqs),
                        exseqs_desc,
                        seq))
                return (None, False)

        liner.send(
            '{}: {:,} {}\n'.format(
                exseqs_field,
                len(exseqs),
                exseqs_desc,))
        df_valid = True
    else:
        # Should not reach here if logic is correct
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                exseqs_field,
                exseqs))
        return (None, False)

    # Return Results
    return (exseqs, df_valid)

def _get_exmotif_source_label(source, idx):
    '''
    Generate a display label for an exmotif source.
    Internal use only.

    :: source
       type - string / list / pd.DataFrame / etc.
       desc - the raw source value
    :: idx
       type - int
       desc - 1-based index for synthetic labels
    '''

    if isinstance(source, str):
        if ',' in source:
            return 'comma[{}]'.format(idx)
        if os.sep in source or '.' in source:
            return _display_path(source)
        return 'inline[{}]'.format(idx)
    elif isinstance(source, pd.DataFrame):
        return 'DataFrame[{}]'.format(idx)
    elif isinstance(source, (list, tuple, set)):
        return 'inline[{}]'.format(idx)
    else:
        return 'input[{}]'.format(idx)

def get_parsed_exmotifs(
    exmotifs,
    exmotifs_field,
    liner):
    '''
    Validate single or multiple excluded-motif inputs.
    Returns (exmotifs, exmotifs_valid, exmotif_inputs).
    Internal use only.

    Polymorphic wrapper around get_parsed_exseqs_info.
    Accepts:
      - None                  -> no screening
      - single set            -> list of motifs, comma-string,
                                CSV/FASTA path, DataFrame
      - multiple sets (list)  -> list of set sources
      - named sets  (dict)    -> {name: set_source, ...}

    Disambiguation for flat list[str]:
      - All elements pass strict ATGC -> single motif set
      - Otherwise -> each element is a separate set source

    Returns:
      - exmotifs       : merged union list (or None)
      - exmotifs_valid : boolean
      - exmotif_inputs : list of dicts with per-input
                         breakdown (or None)

    :: exmotifs
       type - list / str / dict / pd.DataFrame / None
       desc - excluded motif input(s)
    :: exmotifs_field
       type - string
       desc - fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    exmotifs_field = _normalize_field(exmotifs_field)

    # None -> no screening
    if exmotifs is None:
        liner.send(
            '{}: 0 Unique Motif(s)\n'.format(
                exmotifs_field))
        return (None, True, None)

    # Pre-process: comma-string splitting (Python API)
    if isinstance(exmotifs, str) and ',' in exmotifs:
        exmotifs = [s.strip()
            for s in exmotifs.split(',') if s.strip()]

    # Classify input shape

    is_multi = False
    sources = []  # list of (name, raw_source)

    if isinstance(exmotifs, dict):
        # Named multi-set
        is_multi = True
        sources = [(k, v) for k, v in exmotifs.items()]

    elif isinstance(exmotifs, (str, pd.DataFrame)):
        # Single set: string (path/DNA) or DataFrame
        is_multi = False

    elif hasattr(exmotifs, '__iter__'):
        # List/iterable - disambiguate
        items = list(exmotifs)

        if not items:
            liner.send(
                '{}: 0 Unique Motif(s)\n'.format(
                    exmotifs_field))
            return (None, True, None)

        # All strict ATGC strings -> single motif set
        if all(ut.is_strict_DNA(s) for s in items):
            is_multi = False
        else:
            is_multi = True
            sources = [(None, item) for item in items]
    else:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                exmotifs_field, exmotifs))
        return (None, False, None)

    # Single-set path

    if not is_multi:
        (motifs, valid) = get_parsed_exseqs_info(
            exseqs=exmotifs,
            exseqs_field=exmotifs_field,
            exseqs_desc='Unique Motif(s)',
            df_field='Exmotif',
            required=False,
            liner=liner)

        if not valid:
            return (None, False, None)
        if motifs is None:
            return (None, True, None)

        # Strict ATGC enforcement (post-validation)
        for m in motifs:
            if not ut.is_strict_DNA(m):
                liner.send(
                    '{}: [NON-ATGC MOTIF=\'{}\']\n'.format(
                        exmotifs_field, m))
                return (None, False, None)

        # Build single exmotif_inputs entry
        if isinstance(exmotifs, str):
            label = _get_exmotif_source_label(
                source=exmotifs, idx=1)
        elif isinstance(exmotifs, pd.DataFrame):
            label = 'DataFrame'
        else:
            label = 'inline'

        exmotif_inputs = [{
            'name':        None,
            'label':       label,
            'motif_count': len(motifs),
            'motifs':      list(motifs),
        }]

        return (motifs, True, exmotif_inputs)

    # Multi-set path

    # Single-item list -> delegate to single-set display
    if len(sources) == 1:
        name, source = sources[0]

        # Pre-process comma-strings
        if isinstance(source, str) and ',' in source:
            source = [s.strip()
                for s in source.split(',') if s.strip()]

        (motifs, valid) = get_parsed_exseqs_info(
            exseqs=source,
            exseqs_field=exmotifs_field,
            exseqs_desc='Unique Motif(s)',
            df_field='Exmotif',
            required=False,
            liner=liner)

        if not valid:
            return (None, False, None)
        if motifs is None:
            return (None, True, None)

        # Strict ATGC enforcement
        for m in motifs:
            if not ut.is_strict_DNA(m):
                liner.send(
                    '{}: [NON-ATGC MOTIF=\'{}\']\n'.format(
                        exmotifs_field, m))
                return (None, False, None)

        label = name if name else _get_exmotif_source_label(
            source=sources[0][1], idx=1)
        exmotif_inputs = [{
            'name':        name,
            'label':       label,
            'motif_count': len(motifs),
            'motifs':      list(motifs),
        }]

        return (motifs, True, exmotif_inputs)

    # Multiple sources - custom display
    altspacing = ' ' * len(exmotifs_field)

    liner.send(
        '{}: {:,} Exmotif Input(s)\n'.format(
            exmotifs_field, len(sources)))

    all_valid      = True
    all_motifs     = []
    exmotif_inputs = []
    seen           = set()

    # Silent liner for per-source validation
    silent = ut.liner_engine(online=False)

    for idx, (name, source) in enumerate(sources, 1):

        # Determine label BEFORE pre-processing
        label = name if name else _get_exmotif_source_label(
            source=source, idx=idx)

        # Pre-process comma-strings within each source
        if isinstance(source, str) and ',' in source:
            source = [s.strip()
                for s in source.split(',') if s.strip()]

        # Duplicate detection (file paths / DataFrame refs)
        dup_key = None
        if isinstance(source, str) and \
           not ut.is_strict_DNA(source):
            dup_key = os.path.abspath(source)
        elif isinstance(source, pd.DataFrame):
            dup_key = id(source)

        if dup_key is not None and dup_key in seen:
            liner.send(
                '{}:  [{}] {} [DUPLICATE EXMOTIF INPUT]\n'.format(
                    altspacing, idx, label))
            all_valid = False
            continue

        if dup_key is not None:
            seen.add(dup_key)

        # Parse this source silently
        (motifs, valid) = get_parsed_exseqs_info(
            exseqs=source,
            exseqs_field=exmotifs_field,
            exseqs_desc='Unique Motif(s)',
            df_field='Exmotif',
            required=False,
            liner=silent)

        if not valid:
            liner.send(
                '{}:  [{}] {} [INVALID EXMOTIF INPUT]\n'.format(
                    altspacing, idx, label))
            all_valid = False
            continue

        if motifs is None or len(motifs) == 0:
            liner.send(
                '{}:  [{}] {} (0 Unique Motif(s))\n'.format(
                    altspacing, idx, label))
            exmotif_inputs.append({
                'name':        name,
                'label':       label,
                'motif_count': 0,
                'motifs':      [],
            })
            continue

        # Per-source strict ATGC enforcement
        source_invalid = False
        for m in motifs:
            if not ut.is_strict_DNA(m):
                liner.send(
                    '{}:  [{}] {} [NON-ATGC MOTIF=\'{}\']\n'.format(
                        altspacing, idx, label, m))
                all_valid = False
                source_invalid = True
                break

        if source_invalid:
            continue

        liner.send(
            '{}:  [{}] {} ({:,} Unique Motif(s))\n'.format(
                altspacing, idx, label, len(motifs)))

        all_motifs.extend(motifs)
        exmotif_inputs.append({
            'name':        name,
            'label':       label,
            'motif_count': len(motifs),
            'motifs':      list(motifs),
        })

    silent.close()

    if not all_valid:
        return (None, False, None)

    # Merge and deduplicate across all sources
    merged = list(ut.get_uniques(
        iterable=all_motifs,
        typer=list))

    # Print merged total
    liner.send(
        '{}:  {:,} Unique Motif(s)\n'.format(
            altspacing, len(merged)))

    if not merged:
        return (None, True, None)

    return (merged, True, exmotif_inputs)

def get_numeric_validity(
    numeric,
    numeric_field,
    numeric_pre_desc,
    numeric_post_desc,
    minval,
    maxval,
    precheck,
    liner):
    '''
    Determine if numeric is a Real number.
    Internal use only.

    :: numeric
       type - Real
       desc - number to validate
    :: numeric_field
       type - string
       desc - numeric fieldname used in
              printing
    :: numeric_pre_desc
       type - string
       desc - numeric pre-description
              used in printing
    :: numeric_post_desc
       type - string
       desc - numeric post-description
              used in printing
    :: minval
       type - Real
       desc - minimum allowed value
    :: maxval
       type - Real
       desc - maximum allowed value
    :: precheck
       type - boolean
       desc - if False prints numeric_desc
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    numeric_field = _normalize_field(numeric_field)

    # numeric is real?
    numeric_valid = False
    if not isinstance(numeric, nu.Real):
        liner.send('{}:{}{}{} [INPUT TYPE IS INVALID]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc))

    # numeric less than lowerbound?
    elif numeric < minval:
        liner.send('{}:{}{:,}{} [INPUT VALUE SMALLER THAN {:,}]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc, minval))

    # numeric greater than upperbound?
    elif numeric > maxval:
        liner.send('{}:{}{:,}{} [INPUT VALUE LARGER THAN {:,}]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc, maxval))

    # All conditions met!
    else:
        if not precheck:
            liner.send('{}:{}{:,}{}\n'.format(
                numeric_field, numeric_pre_desc, numeric, numeric_post_desc))
        numeric_valid = True

    # Return numeric validity
    return numeric_valid

def get_optional_numeric_validity(
    numeric,
    numeric_field,
    numeric_pre_desc,
    numeric_post_desc,
    minval,
    maxval,
    precheck,
    liner):
    '''
    Determine if optionally provided
    numeric is a Real number.
    Internal use only.

    :: numeric
       type - Real
       desc - number to validate
    :: numeric_field
       type - string
       desc - numeric fieldname used in
              printing
    :: numeric_pre_desc
       type - string
       desc - numeric pre-description
              used in printing
    :: numeric_post_desc
       type - string
       desc - numeric post-description
              used in printing
    :: minval
       type - Real
       desc - minimum allowed value
    :: maxval
       type - Real
       desc - maximum allowed value
    :: precheck
       type - boolean
       desc - if False prints numeric_desc
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # numeric is None
    if numeric is None:
        liner.send(
            '{}: None Specified\n'.format(
                numeric_field))
        return True

    # Regular numeric Validation
    return get_numeric_validity(
        numeric=numeric,
        numeric_field=numeric_field,
        numeric_pre_desc=numeric_pre_desc,
        numeric_post_desc=numeric_post_desc,
        minval=minval,
        maxval=maxval,
        precheck=precheck,
        liner=liner)

def get_sample_size_validity(
    sample_size,
    sample_size_field,
    liner):
    '''
    Determine if sample_size is a valid positive integer.
    Internal use only.

    :: sample_size
       type - integer
       desc - sample size to validate
    :: sample_size_field
       type - string
       desc - sample size fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    return get_numeric_validity(
        numeric=sample_size,
        numeric_field=sample_size_field,
        numeric_pre_desc=' ',
        numeric_post_desc=' Samples per Category',
        minval=1,
        maxval=100000,
        precheck=False,
        liner=liner)

def get_parsed_spacerlen_info(
    spacerlen,
    spacerlen_field,
    df_field,
    oligolimit,
    oligolimit_valid,
    indf,
    indata_valid,
    liner):
    '''
    Determine if given spacer length(s) input are valid.
    Internal use only.

    :: spacerlen
       type - iterable / integer / pd.DataFrame / None
       desc - iterable of integers denoting the length
              of spacers; optionally a DataFrame or a
              path to a CSV file with such lengths
    :: spacerlen_field
       type - string
       desc - spacerlen fieldname used in printing
    :: df_field
       type - string
       desc - name of the column in DataFrame storing
              spacer lengths
    :: oligolimit
       type - integer
       desc - maximum oligo length after inserting
              designed spacers
    :: oligolimit_valid
       type - boolean
       desc - oligolimit parsing status
    :: indf
       type - pd.DataFrame / None
       desc - Associated input DataFrame
    :: indata_valid
       type - boolean
       desc - Associated input DataFrame validity
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    spacerlen_field = _normalize_field(spacerlen_field)
    df_field = _normalize_field(df_field)

    # Is spacerlen None?
    if spacerlen is None:
        liner.send(
            '{}: Computed from Oligo Length (Auto-Inferred)\n'.format(
                spacerlen_field))
        return (None, True)

    # Quick compute indexlen
    if indata_valid:
        lenindex = len(indf.index)
    else:
        lenindex = 1

    # Correct oligolimit
    if not oligolimit_valid:
        oligolimit = float('+inf')
    else:
        oligolimit = round(oligolimit)

    # Is spacerlen numeric?
    if isinstance(spacerlen, nu.Real):

        # Ensure spacerlen valid
        spacerlen_valid = get_numeric_validity(
            numeric=spacerlen,
            numeric_field=spacerlen_field,
            numeric_pre_desc=' Exactly ',
            numeric_post_desc=' Base Pair(s)',
            minval=1,
            maxval=oligolimit,
            precheck=False,
            liner=liner)

        # Build spacerlen
        if spacerlen_valid:
            spacerlen = np.zeros(
                lenindex, dtype=np.int64) + round(spacerlen)
        else:
            spacerlen = None

        # Return Results
        return (spacerlen, spacerlen_valid)

    # Is spacerlen extracted?
    spacerlen_extracted = False

    # Is spacerlen iterable?
    if not isinstance(spacerlen, pd.DataFrame) and \
       not isinstance(spacerlen, str) and \
       not isinstance(spacerlen, nu.Real):

        # Try extracting spacerlen
        try:
            spacerlen = list(n for n in spacerlen)

        # Error during extraction
        except:
            # Handled in sink below:
            # 'if not spacerlen_extracted: ...'
            spacerlen_extracted = False

        # Iteration successful
        else:
            spacerlen_extracted = True

        # Do we have enough spacers?
        if indata_valid:
            if len(spacerlen) != len(indf.index):

                # Compute category of mismatch
                catm = ['FEWER', 'MORE'][len(spacerlen) > len(indf.index)]

                # Show Update
                liner.send(
                    '{}: Iterable of {:,} Record(s) [{} THAN {:,} SPACERS SPECIFIED]\n'.format(
                        spacerlen_field,
                        len(spacerlen),
                        catm,
                        len(indf.index)))

                # Return Results
                return (None, False)

    # Is spacerlen a CSV file or DataFrame?
    else:

        # Compute spacerlen data validity
        (df,
        data_name,
        df_valid) = get_parsed_data_info(
            data=spacerlen,
            data_field=spacerlen_field,
            required_fields=('ID', df_field,),
            liner=liner)

        # Is spacerlen df valid?
        if df_valid:
            # Are the indexes matching?
            idx_match = False
            if indata_valid:

                # Matching IDs?
                idx_match = len(df.index)    == len(indf.index) and \
                            sorted(df.index) == sorted(indf.index)

                # Everything OK!
                if idx_match:
                    df = df.reindex(indf.index)
                    spacerlen = list(n for n in df[df_field])
                    spacerlen_extracted = True

            # No basis for matching, so show indifference
            # here, return and fail later on
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s)\n'.format(
                        spacerlen_field,
                        data_name,
                        len(df.index)))
                return (None, False)

            # Indexes don't match
            if not idx_match:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH INPUT DATA]\n'.format(
                        spacerlen_field,
                        data_name,
                        len(df.index)))
                return (None, False)

        else:
            # Note: We've already shown how
            # invalidity impacts spacerlen,
            # so we're directly returning
            return (None, False)

    # spacerlen not extracted?
    if not spacerlen_extracted:
        # Note: if spacerlen is neither numeric, pd.DataFrame
        # iterable, nor a string path, then it should sink here
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                spacerlen_field,
                spacerlen))
        return (None, False)

    # Finalize spacerlen
    minimum = float('+inf')
    maximum = float('-inf')
    finalspacerlen = []

    for sl in spacerlen:

        # Ensure sl valid
        sl_valid = get_numeric_validity(
            numeric=sl,
            numeric_field=spacerlen_field,
            numeric_pre_desc=' Found an Entry for ',
            numeric_post_desc=' Base Pair(s)',
            minval=0,
            maxval=oligolimit,
            precheck=True,
            liner=liner)

        # Invalid sl
        if not sl_valid:
            return (None, False)
        # So far valid
        else:
            sl = round(sl)
            finalspacerlen.append(sl)
            minimum = min(minimum, sl)
            maximum = max(maximum, sl)

    # No errors whatsoever
    spacerlen = finalspacerlen
    spacerlen_valid = True

    if minimum == maximum:

        # We don't tolerate all zeros
        if minimum == 0:
            liner.send(
                '{}: Exactly 0 Base Pair(s) [ALL SPACERS HAVE ZERO LENGTH]\n'.format(
                    spacerlen_field,
                    minimum))
            return (None, False)

        # At least some entries require spacers
        else:
            liner.send(
                '{}: Exactly {:,} Base Pair(s)\n'.format(
                    spacerlen_field,
                    minimum))
    else:
        liner.send(
            '{}: {:,} to {:,} Base Pair(s)\n'.format(
                spacerlen_field,
                minimum,
                maximum))

    # Pack spacerlen
    spacerlen = np.array(spacerlen)

    # Return Results
    return (spacerlen, True)

def get_parsed_oligosets_info(
    oligosets,
    oligosets_field,
    df_field,
    indf,
    indata_valid,
    liner):
    '''
    Determine if oligosets input is valid.
    Internal use only.

    :: oligosets
       type - iterable / string / pd.DataFrame / None
       desc - iterable of set labels per oligo; optionally
              a DataFrame or a path to a CSV file with such
              labels
    :: oligosets_field
       type - string
       desc - oligosets fieldname used in printing
    :: df_field
       type - string
       desc - name of the column in DataFrame storing
              oligo set labels
    :: indf
       type - pd.DataFrame / None
       desc - associated input DataFrame
    :: indata_valid
       type - boolean
       desc - associated input DataFrame validity
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    oligosets_field = _normalize_field(oligosets_field)

    # Is oligosets None?
    if oligosets is None:
        liner.send(
            '{}: None Specified\n'.format(
                oligosets_field))
        return (None, True)

    # Quick compute indexlen
    if indata_valid:
        lenindex = len(indf.index)
    else:
        lenindex = 1

    # Is oligosets iterable?
    if not isinstance(oligosets, pd.DataFrame) and \
       not isinstance(oligosets, str):

        # Try extracting oligosets
        try:
            oligosets = list(o for o in oligosets)

        # Error during extraction
        except:
            liner.send(
                '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                    oligosets_field,
                    oligosets))
            return (None, False)

        # Do we have enough set labels?
        if indata_valid:
            if len(oligosets) != len(indf.index):

                # Compute category of mismatch
                catm = ['FEWER', 'MORE'][len(oligosets) > len(indf.index)]

                # Show Update
                liner.send(
                    '{}: Iterable of {:,} Record(s) [{} THAN {:,} SET LABELS SPECIFIED]\n'.format(
                        oligosets_field,
                        len(oligosets),
                        catm,
                        len(indf.index)))

                # Return Results
                return (None, False)

        # No basis for matching, so show indifference
        # here, return and fail later on
        else:
            liner.send(
                '{}: Iterable of {:,} Record(s)\n'.format(
                    oligosets_field,
                    len(oligosets)))
            return (None, False)

    # Is oligosets a CSV file or DataFrame?
    else:

        # Compute oligosets data validity
        (df,
        data_name,
        df_valid) = get_parsed_data_info(
            data=oligosets,
            data_field=oligosets_field,
            required_fields=('ID', df_field,),
            liner=liner)

        # Is oligosets df valid?
        if df_valid:

            # Are the indexes matching?
            idx_match = False
            if indata_valid:

                # Matching IDs?
                idx_match = len(df.index)    == len(indf.index) and \
                            sorted(df.index) == sorted(indf.index)

                # Everything OK!
                if idx_match:
                    df = df.reindex(indf.index)
                    oligosets = list(o for o in df[df_field])

            # No basis for matching, so show indifference
            # here, return and fail later on
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s)\n'.format(
                        oligosets_field,
                        data_name,
                        len(df.index)))
                return (None, False)

            # Indexes don't match
            if not idx_match:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH INPUT DATA]\n'.format(
                        oligosets_field,
                        data_name,
                        len(df.index)))
                return (None, False)

        else:
            # Note: We've already shown how
            # invalidity impacts oligosets,
            # so we're directly returning
            return (None, False)

    # Validate set labels (hashable, non-missing)
    for label in oligosets:

        # Unhashable label?
        try:
            hash(label)
        except TypeError:
            liner.send(
                '{}: {:,} Record(s) [UNHASHABLE SET LABEL=\'{}\']\n'.format(
                    oligosets_field,
                    len(oligosets),
                    label))
            return (None, False)

        # Missing label?
        if pd.isna(label):
            liner.send(
                '{}: {:,} Record(s) [MISSING SET LABEL]\n'.format(
                    oligosets_field,
                    len(oligosets)))
            return (None, False)

    # Count unique sets
    uniques = ut.get_uniques(
        iterable=oligosets,
        typer=tuple)
    liner.send(
        '{}: {:,} Unique Set(s)\n'.format(
            oligosets_field,
            len(uniques)))

    # Return Results
    return (np.array(oligosets, dtype=object), True)

def get_parsed_pairedprimer_info(
    pairedcol,
    pairedcol_field,
    df,
    oligosets,
    allow_missing,
    liner):
    '''
    Determine if paired primer input is valid,
    including per-set pairing when oligosets
    is provided. Internal use only.

    :: pairedcol
       type - string / None
       desc - column name in df storing paired
              primer sequences
    :: pairedcol_field
       type - string
       desc - pairedcol fieldname used in printing
    :: df
       type - pd.DataFrame / None
       desc - input DataFrame where paired primers
              are present
    :: oligosets
       type - np.array / None
       desc - per-row oligo set labels
    :: allow_missing
       type - boolean
       desc - if True, allow sets with all
              missing paired primers
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    pairedcol_field = _normalize_field(pairedcol_field)

    # Is pairedcol None?
    if pairedcol is None:
        liner.send(
            '{}: None Specified\n'.format(
                pairedcol_field))
        return None, True

    # Is pairedcol a string?
    if not isinstance(pairedcol, str):
        liner.send(
            '{}: Input from Column \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                pairedcol_field,
                pairedcol))
        return None, False

    # Do we have a DataFrame?
    if not isinstance(df, pd.DataFrame):
        liner.send(
            '{}: Input from Column \'{}\' [COLUMN EXISTENCE UNKNOWN]\n'.format(
                pairedcol_field,
                pairedcol))
        return None, False

    # Column missing?
    if pairedcol not in df.columns:
        liner.send(
            '{}: Input from Column \'{}\' [MISSING COLUMN=\'{}\']\n'.format(
                pairedcol_field,
                pairedcol,
                pairedcol))
        return None, False

    # No oligosets: require a single constant primer
    if oligosets is None and not allow_missing:
        return get_constantcol_validity(
            constantcol=pairedcol,
            constantcol_field=pairedcol_field,
            df=df,
            liner=liner)

    # Missing-aware parsing
    if allow_missing:

        # No oligosets: allow all-missing column
        if oligosets is None:
            missing_mask = ut.get_missing_mask(
                series=df[pairedcol],
                allow_dash=True)
            if missing_mask.all():
                liner.send(
                    '{}: None Specified\n'.format(
                        pairedcol_field))
                return None, True
            if missing_mask.any():
                liner.send(
                    '{}: Input from Column \'{}\' [MISSING VALUE(S)]\n'.format(
                        pairedcol_field,
                        pairedcol))
                return None, False

            # Extract unique candidates
            uniques = ut.get_uniques(
                iterable=df[pairedcol],
                typer=tuple)
            if len(uniques) > 1:
                liner.send(
                    '{}: Input from Column \'{}\' [NON-UNIQUE COLUMN=\'{}\']\n'.format(
                        pairedcol_field,
                        pairedcol,
                        pairedcol))
                return None, False

            liner.send(
                '{}: A {:,} Base Pair DNA Sequence\n'.format(
                    pairedcol_field,
                    len(uniques[0])))
            return uniques[0], True

    # Oligosets provided: require constant per set
    (uniques,
    groups,
    _) = ut.get_oligoset_groups(
        oligosets=oligosets)
    paired_map = {}

    for label in uniques:
        idx = groups[label]
        values = df[pairedcol].iloc[idx]
        missing_mask = ut.get_missing_mask(
            series=values,
            allow_dash=True)
        values = values[~missing_mask]

        # Missing-only set?
        if allow_missing and values.empty:
            paired_map[label] = None
            continue

        # Missing with values?
        if allow_missing and missing_mask.any():
            liner.send(
                '{}: Input from Column \'{}\' [MISSING VALUE IN SET=\'{}\']\n'.format(
                    pairedcol_field,
                    pairedcol,
                    label))
            return None, False

        values = ut.get_uniques(
            iterable=values,
            typer=tuple)

        # Non-unique per set?
        if len(values) != 1:
            liner.send(
                '{}: Input from Column \'{}\' [NON-UNIQUE WITHIN SET=\'{}\']\n'.format(
                    pairedcol_field,
                    pairedcol,
                    label))
            return None, False

        paired_map[label] = values[0]

    liner.send(
        '{}: Input from Column \'{}\' (Per-Set)\n'.format(
            pairedcol_field,
            pairedcol))

    return paired_map, True

def get_parsed_compress_data_info(
    data,
    data_field,
    liner):
    '''
    Parse input for compress: validate ID column exists and all other columns
    are non-empty DNA sequences (A/T/G/C only, no degenerate codes).
    Internal use only.

    :: data
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas DataFrame storing variant sequences
    :: data_field
       type - string
       desc - data fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing

    Returns:
        Tuple of (df, dna_columns, valid):
        - df: parsed DataFrame indexed by ID (or None if invalid)
        - dna_columns: list of DNA sequence column names
        - valid: boolean indicating overall validity
    '''

    data_field = _normalize_field(data_field)

    # Is data valid CSV or DataFrame?
    (df,
    data_name,
    df_valid) = get_parsed_data_info(
        data=data,
        data_field=data_field,
        required_fields=('ID',),
        liner=liner)

    dna_columns = []

    if not df_valid:
        return (None, dna_columns, False)

    # Get all non-ID columns as potential DNA columns
    dna_columns = [col for col in df.columns]

    if not dna_columns:
        liner.send(
            '{}: {} w/ {:,} Record(s) [NO DNA COLUMNS FOUND]\n'.format(
                data_field,
                data_name,
                len(df.index)))
        return (None, [], False)

    # Check that all entries are concrete DNA (A/T/G/C only)
    strict_dna_alpha = set('ATGC')

    for column in dna_columns:
        for idx, value in enumerate(df[column]):
            if not isinstance(value, str):
                liner.send(
                    '{}: {} w/ {:,} Record(s) [NON-STRING VALUE IN COLUMN=\'{}\']\n'.format(
                        data_field,
                        data_name,
                        len(df.index),
                        column))
                return (None, [], False)

            value_upper = value.upper().strip()
            if not value_upper:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [EMPTY VALUE IN COLUMN=\'{}\']\n'.format(
                        data_field,
                        data_name,
                        len(df.index),
                        column))
                return (None, [], False)

            if not set(value_upper) <= strict_dna_alpha:
                invalid_chars = set(value_upper) - strict_dna_alpha
                examples = ut.get_row_examples(
                    df=df,
                    invalid_mask=~df[column].str.upper().str.match(r'^[ATGC]+$'),
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '{}: {} w/ {:,} Record(s) [NON-CONCRETE DNA (A/T/G/C only) IN COLUMN=\'{}\']{}\n'.format(
                        data_field,
                        data_name,
                        len(df.index),
                        column,
                        example_note))
                return (None, [], False)

    # Uppercase all DNA columns
    for col in dna_columns:
        df[col] = df[col].str.upper()

    # Show success
    liner.send(
        '{}: {} w/ {:,} Variant(s)\n'.format(
            data_field,
            data_name,
            len(df.index)))

    return (df, dna_columns, True)

def get_parsed_expand_data_info(
    data,
    data_field,
    sequence_column,
    sequence_column_field,
    liner):
    '''
    Parse input for expand: validate ID + sequence_column exist and
    sequence_column contains valid IUPAC sequences.
    Internal use only.

    :: data
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas DataFrame storing degenerate sequences
    :: data_field
       type - string
       desc - data fieldname used in printing
    :: sequence_column
       type - string
       desc - column name containing IUPAC sequences
    :: sequence_column_field
       type - string
       desc - sequence_column fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing

    Returns:
        Tuple of (df, valid):
        - df: parsed DataFrame indexed by ID (or None if invalid)
        - valid: boolean indicating overall validity
    '''

    data_field = _normalize_field(data_field)
    sequence_column_field = _normalize_field(sequence_column_field)

    # Is data valid CSV or DataFrame?
    (df,
    data_name,
    df_valid) = get_parsed_data_info(
        data=data,
        data_field=data_field,
        required_fields=('ID',),
        liner=liner,
        id_aliases=('DegenerateID',))

    if not df_valid:
        return (None, False)

    # Show input data success
    liner.send(
        '{}: {} w/ {:,} Oligo(s)\n'.format(
            data_field,
            data_name,
            len(df.index)))

    # Check sequence_column is a valid string
    if not isinstance(sequence_column, str):
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                sequence_column_field, sequence_column))
        return (None, False)

    # Check sequence_column exists
    if sequence_column not in df.columns:
        liner.send(
            '{}: \'{}\' [COLUMN DOES NOT EXIST]\n'.format(
                sequence_column_field, sequence_column))
        return (None, False)

    # Check that all values in sequence_column are valid IUPAC sequences
    iupac_alpha = set(ut.ddna_space.keys()) - {'-'}

    for idx, value in enumerate(df[sequence_column]):
        if not isinstance(value, str):
            liner.send(
                '{}: \'{}\' [NON-STRING VALUE IN COLUMN]\n'.format(
                    sequence_column_field, sequence_column))
            return (None, False)

        value_upper = value.upper().strip()
        if not value_upper:
            liner.send(
                '{}: \'{}\' [EMPTY VALUE IN COLUMN]\n'.format(
                    sequence_column_field, sequence_column))
            return (None, False)

        if not set(value_upper) <= iupac_alpha:
            invalid_chars = set(value_upper) - iupac_alpha
            liner.send(
                '{}: \'{}\' [INVALID IUPAC CHARACTER(S): {}]\n'.format(
                    sequence_column_field, sequence_column, sorted(invalid_chars)))
            return (None, False)

    # Uppercase the sequence column
    df[sequence_column] = df[sequence_column].str.upper()

    # Show success
    liner.send(
        '{}: \'{}\' ... Parsed\n'.format(
            sequence_column_field, sequence_column))

    return (df, True)

def get_parsed_compress_mapping_info(
    mapping_data,
    mapping_field,
    liner,
    required=False):
    '''
    Parse mapping output from `compress` for use in `expand`.
    Supports a CSV path (auto-suffixed) or a DataFrame. Unlike `get_parsed_data_info`, preserves `ID`
    as a normal column (not as the index), since downstream joins expect it.
    Internal use only.

    :: mapping_data
       type - string / pd.DataFrame / None
       desc - path to mapping CSV output by `compress` or a mapping DataFrame
              (expected suffix='.oligopool.compress.mapping.csv')
    :: mapping_field
       type - string
       desc - mapping fieldname used in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    :: required
       type - bool
       desc - if True, mapping_data must be provided (default: False)
    '''

    mapping_field = _normalize_field(mapping_field)

    # None specified
    if mapping_data is None:
        if required:
            liner.send('{}: None Specified [REQUIRED]\n'.format(mapping_field))
            return (None, False)
        liner.send('{}: None Specified\n'.format(mapping_field))
        return (None, True)

    # DataFrame specified
    if isinstance(mapping_data, pd.DataFrame):
        mapdf = mapping_data.copy()
        # If ID is provided as the index (common in library usage), restore it as a column.
        if 'ID' not in mapdf.columns:
            if mapdf.index.name == 'ID':
                mapdf = mapdf.reset_index()
            else:
                mapdf = mapdf.reset_index(drop=False)
                if 'index' in mapdf.columns and 'ID' not in mapdf.columns:
                    mapdf.rename(columns={'index': 'ID'}, inplace=True)
        liner.send(
            '{}: DataFrame w/ {:,} Row(s)\n'.format(
                mapping_field,
                len(mapdf.index)))
    # Path specified
    elif isinstance(mapping_data, str):
        # Validate + auto-suffix
        mapping_ok = get_infile_validity(
            infile=mapping_data,
            infile_suffix='.oligopool.compress.mapping.csv',
            infile_field=mapping_field,
            liner=liner)
        if not mapping_ok:
            return (None, False)

        mapfile = ut.get_adjusted_path(
            path=mapping_data,
            suffix='.oligopool.compress.mapping.csv')
        try:
            mapdf = pd.read_csv(mapfile, sep=',', header=0, engine='c')
        except Exception:
            liner.send('{}: {} [INVALID CSV FILE]\n'.format(
                mapping_field,
                _display_path(mapfile)))
            return (None, False)
        liner.send('{}: Loaded from {}\n'.format(
            mapping_field,
            _display_path(mapfile)))
    else:
        liner.send('{}: {} [INPUT TYPE IS INVALID]\n'.format(mapping_field, mapping_data))
        return (None, False)

    # Required columns
    required_cols = ('ID', 'DegenerateID')
    missing_cols = [c for c in required_cols if c not in mapdf.columns]
    if missing_cols:
        liner.send(
            '{}: [MISSING COLUMN(S)={}] \n'.format(
                mapping_field,
                ','.join(missing_cols)))
        return (None, False)

    # Required columns must be non-missing
    for col in required_cols:
        missing_mask = ut.get_missing_mask(series=mapdf[col], allow_dash=False)
        if missing_mask.any():
            examples = ut.get_row_examples(
                df=mapdf,
                invalid_mask=missing_mask,
                id_col='ID',
                limit=5)
            example_note = ut.format_row_examples(examples)
            liner.send(
                '{}: [MISSING VALUES IN COLUMN=\'{}\']{}\n'.format(
                    mapping_field,
                    col,
                    example_note))
            return (None, False)

    # Normalize key columns to strings for consistent merges.
    mapdf['ID'] = mapdf['ID'].map(str)
    mapdf['DegenerateID'] = mapdf['DegenerateID'].map(str)

    return (mapdf, True)

def get_parsed_associatedata_info(
    associatedata,
    associatedata_field,
    required_fields,
    bardf,
    barcodedata_valid,
    liner):
    '''
    Determine if associatedata is valid and the
    ID is consistant with bardf.
    Internal use only.

    :: associatedata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing variant
              information
    :: associatedata_field
       type - string
       desc - associatedata fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: bardf
       type - pd.DataFrame / None
       desc - pandas DataFrame containing
              barcode information
    :: barcodedata_valid
       type - boolean
       desc - if True indicates that shared
              ID is to be computed, otherwise
              ID matching is skipped
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is associatedata None?
    if associatedata is None:
        liner.send(
            '{}: None Specified\n'.format(
                associatedata_field))
        return (None, True) # Because variant is optional

    # Is associatedata valid?
    (assdf,
    data_name,
    associatedata_valid) = get_parsed_indata_info(
        indata=associatedata,
        indata_field=associatedata_field,
        required_fields=required_fields,
        precheck=True,
        liner=liner)

    # Does assdf share ID with bardf?
    associatedata_idx_match = False
    if associatedata_valid:

        # Are the indexes matching?
        idx_match = False
        if barcodedata_valid:

            # Matching IDs?
            idx_match = len(assdf.index)    == len(bardf.index) and \
                        sorted(assdf.index) == sorted(bardf.index)

            # Everything OK!
            if idx_match:
                assdf = assdf.reindex(bardf.index)
                associatedata_idx_match = True # same as idx_match

            # Indexes don't match
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH BARCODE DATA]\n'.format(
                        associatedata_field,
                        data_name,
                        len(assdf.index)))

        # No basis for matching, so show indifference
        # here, return and fail later on
        else:
            associatedata_idx_match = True # Technically, we can't complain!

    # Compute final validity
    assdf_valid = all([
        associatedata_valid,
        associatedata_idx_match])

    # Is df valid?
    if assdf_valid:
        liner.send(
            '{}: {} w/ {:,} Record(s)\n'.format(
                associatedata_field,
                data_name,
                len(assdf.index)))
    else:
        # Erase df
        assdf = None

    # Return data validity
    return (assdf, assdf_valid)

def resolve_type_parameter(
    value,
    type_name):
    '''
    Resolve a type parameter value (int or string)
    to its canonical integer value using the
    TYPE_REGISTRIES. Internal use only.

    :: value
       type - int or string
       desc - the type parameter value to
              resolve
    :: type_name
       type - string
       desc - registry key (e.g. 'barcode_type',
              'primer_type', 'motif_type')

    Returns (resolved_int, valid_bool, warning_msg)
    tuple where resolved_int is the canonical integer,
    valid_bool indicates success, and warning_msg
    contains any fuzzy match warning or None.
    '''

    registry = TYPE_REGISTRIES.get(type_name)
    if registry is None:
        return (None, False, None)

    aliases = registry['aliases']

    # Integer pass-through (accept numpy integer types; reject bool)
    if isinstance(value, bool):
        return (None, False, None)
    if isinstance(value, nu.Integral):
        ival = int(value)
        if ival in aliases:
            return (aliases[ival], True, None)
        return (None, False, None)

    # String resolution
    if isinstance(value, str):
        normalized = value.strip().lower()

        # Numeric strings ("0", "1") support
        if normalized.lstrip('-').isdigit():
            ival = int(normalized)
            if ival in aliases:
                return (aliases[ival], True, None)

        # Exact match
        if normalized in aliases:
            return (aliases[normalized], True, None)

        # Fuzzy matching
        string_aliases = [k for k in aliases.keys() if isinstance(k, str)]
        matches = difflib.get_close_matches(normalized, string_aliases, n=1, cutoff=0.6)
        if matches:
            matched = matches[0]
            resolved = aliases[matched]
            warning = "Interpreted '{}' as '{}'".format(value, matched)
            return (resolved, True, warning)

        return (None, False, None)

    return (None, False, None)

def get_typed_categorical_validity(
    category,
    category_field,
    category_pre_desc,
    category_post_desc,
    type_name,
    liner):
    '''
    Validate a typed category parameter accepting
    int or string values with alias support and
    fuzzy matching. Internal use only.

    :: category
       type - int or string
       desc - category to validate
    :: category_field
       type - string
       desc - category fieldname used in
              printing
    :: category_pre_desc
       type - string
       desc - category pre-description
              used in printing
    :: category_post_desc
       type - string
       desc - category post-description
              used in printing
    :: type_name
       type - string
       desc - registry key (e.g. 'barcode_type',
              'primer_type', 'motif_type')
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    category_field = _normalize_field(category_field)

    registry = TYPE_REGISTRIES.get(type_name)
    if registry is None:
        liner.send('{}:{}{}{} [UNKNOWN TYPE REGISTRY]\n'.format(
            category_field,
            category_pre_desc,
            category,
            category_post_desc))
        return (category, False)

    # Resolve the value
    resolved, valid, warning = resolve_type_parameter(category, type_name)

    if not valid:
        # Build helpful error message with valid options
        canonical = registry['canonical']
        valid_options = ', '.join(["'{}'".format(c) for c in canonical])
        liner.send('{}:{}{}{} [INVALID: use 0, 1, or {}]\n'.format(
            category_field,
            category_pre_desc,
            category,
            category_post_desc,
            valid_options))
        return (category, False)

    # Show warning if fuzzy matching was used
    if warning:
        liner.send('{}:{}{}{} ({})\n'.format(
            category_field,
            category_pre_desc,
            registry['display'][resolved],
            category_post_desc,
            warning))
    else:
        liner.send('{}:{}{}{}\n'.format(
            category_field,
            category_pre_desc,
            registry['display'][resolved],
            category_post_desc))

    return (resolved, True)

def get_optional_typed_categorical_validity(
    category,
    category_field,
    category_pre_desc,
    category_post_desc,
    type_name,
    liner):
    '''
    Validate an optional typed category parameter
    accepting int or string. If category is None,
    returns (None, True). Internal use only.

    :: category
       type - int, string, or None
       desc - category to validate (None is
              valid for optional params)
    :: category_field
       type - string
       desc - category fieldname used in
              printing
    :: category_pre_desc
       type - string
       desc - category pre-description
              used in printing
    :: category_post_desc
       type - string
       desc - category post-description
              used in printing
    :: type_name
       type - string
       desc - registry key (e.g. 'read_type',
              'pack_type', 'mapping_type')
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    category_field = _normalize_field(category_field)

    # category is None - valid optional
    if category is None:
        liner.send(
            '{}: None Specified\n'.format(
                category_field))
        return (None, True)

    # Regular typed validation
    return get_typed_categorical_validity(
        category=category,
        category_field=category_field,
        category_pre_desc=category_pre_desc,
        category_post_desc=category_post_desc,
        type_name=type_name,
        liner=liner)

def get_parsed_typeIIS_info(
    typeIIS,
    typeIIS_field,
    liner):
    '''
    Determine if typeIIS system selected
    is valid. Internal use only.

    :: typeIIS
       type - string
       desc - name of the TypeIIS enzyme
              selected for excision
    :: typeIIS_field
       type - string
       desc - typeIIS fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is typeIIs string?
    typeIIS_is_string = False
    if not isinstance(typeIIS, str):
        liner.send('{}: Enzyme \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                typeIIS_field,
                typeIIS))
    else:
        typeIIS_is_string = True

    # Is typeIIS known?
    typeIIS_is_known = False
    if typeIIS_is_string:
        typeIIS_ = typeIIS.lower()
        if not typeIIS_ in ut.typeIIS_dict:
            liner.send('{}: Enzyme \'{}\' [UNSUPPORTED ENZYME]\n'.format(
                    typeIIS_field,
                    typeIIS))
        else:
            typeIIS_is_known = True

    # Compute validity
    typeIIS_valid = typeIIS_is_string and typeIIS_is_known

    # Show update
    if typeIIS_valid:
        liner.send('{}: Enzyme \'{}\' Recognizing Motif \'{}\'\n'.format(
            typeIIS_field,
            ut.typeIIS_dict[typeIIS_][0],
            ut.typeIIS_dict[typeIIS_][1]))

    # Return typeIIS validity
    if typeIIS_valid:
        typeIISname = ut.typeIIS_dict[typeIIS_][0]
        typeIIS = ut.typeIIS_dict[typeIIS_][1] + \
                 ('N' * ut.typeIIS_dict[typeIIS_][2])
        return typeIIS, typeIISname, typeIIS_valid
    else:
        return None, None, typeIIS_valid

def get_seqconstr_validity(
    seqconstr,
    seqconstr_field,
    minlenval,
    element,
    liner):
    '''
    Determine if seqconstr is a valid
    degenerate sequence constraint.
    Internal use only.

    :: seqconstr
       type - string
       desc - primer sequence constraint
    :: seqconstr_field
       type - string
       desc - seqconstr fieldname used
              in printing
    :: minlenval
       type - integer
       desc - minimum length allowed for
              sequence constraint
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    seqconstr_field = _normalize_field(seqconstr_field)

    # Is seqconstr string?
    seqconstr_is_string = False
    if not isinstance(seqconstr, str):
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                seqconstr_field,
                seqconstr))
    else:
        seqconstr_is_string = True

    # Is seqconstr degenerate DNA string?
    seqconstr_is_ddna = False
    if seqconstr_is_string:
        if not ut.is_DNA(
            seq=seqconstr,
            dna_alpha=ut.ddna_alpha):
            liner.send(
                '{}: {} [NON-IUPAC VALUE]\n'.format(
                    seqconstr_field,
                    seqconstr))
        else:
            seqconstr_is_ddna = True

    # Is seqconstr long enough?
    seqconstr_is_long = False
    if seqconstr_is_ddna:
        if len(seqconstr) < minlenval:
            liner.send(
                '{}: A {:,} Base Pair IUPAC Constraint [{} SHORTER THAN {:,} BASE PAIRS]\n'.format(
                    seqconstr_field,
                    len(seqconstr),
                    element,
                    minlenval))
        else:
            seqconstr_is_long = True

    # Compute final validity
    seqconstr_valid = all([
        seqconstr_is_string,
        seqconstr_is_ddna,
        seqconstr_is_long])

    # Show update
    if seqconstr_valid:
        liner.send(
            '{}: A {:,} Base Pair IUPAC Constraint\n'.format(
                seqconstr_field,
                len(seqconstr)))

    # Return validity
    return seqconstr_valid

def get_constantcol_validity(
    constantcol,
    constantcol_field,
    df,
    liner):
    '''
    Determine if constantcol is valid.
    Internal use only.

    :: constantcol
       type - string / None
       desc - the name of column in df where
              constantcol is stored
    :: constantcol_field
       type - string
       desc - constantcol fieldname used
              in printing
    :: df
       type - pd.DataFrame
       desc - input DataFrame where constant
              is present
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    constantcol_field = _normalize_field(constantcol_field)

    # Is constantcol None?
    if constantcol is None:
        liner.send(
            '{}: None Specified\n'.format(
                constantcol_field))
        return None, True

    # Is constantcol a string?
    if isinstance(constantcol, str):

        # Is constantcol a column in df?
        if constantcol in df.columns:

            # Extract unique candidates
            uniques = ut.get_uniques(
                iterable=df[constantcol],
                typer=tuple)

            # Too many constant primer candidates?
            if len(uniques) > 1:
                liner.send(
                    '{}: Input from Column \'{}\' [NON-UNIQUE COLUMN=\'{}\']\n'.format(
                        constantcol_field,
                        constantcol,
                        constantcol))
                return None, False

            # Unique constant primer
            else:
                liner.send(
                    '{}: A {:,} Base Pair DNA Sequence\n'.format(
                        constantcol_field,
                        len(uniques[0])))
                return uniques[0], True

        # Nothing matches
        else:
            liner.send(
                '{}: Input from Column \'{}\' [MISSING COLUMN=\'{}\']\n'.format(
                    constantcol_field,
                    constantcol,
                    constantcol))
            return None, False

    # Non-string constantcol
    else:
        liner.send(
            '{}: Input from Column \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                constantcol_field,
                constantcol))
        return None, False

def get_parsed_range_info(
    minval,
    maxval,
    range_field,
    range_unit,
    range_min,
    range_max,
    liner):
    '''
    Determine if range information
    is valid. Internal use only.

    :: minval
       type - nu.Real
       desc - minimum range value
    :: maxval
       type - nu.Real
       desc - maximum range value
    :: range_field
       type - string
       desc - range fieldname used in
              printing
    :: range_unit
       type - string
       desc - range value unit used in
              printing
    :: range_min
       type - nu.Real
       desc - allowed range lowerbound
    :: range_max
       type - nu.Real
       desc - allowed range upperbound
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    range_field = _normalize_field(range_field)

    # Define value range
    vrange = (minval, maxval)

    # Is vrange numeric?
    vrange_is_numeric = False
    if  not isinstance(minval, nu.Real) or \
        not isinstance(maxval, nu.Real):
        liner.send(
            '{}: {} to {} {} [INPUT TYPE IS INVALID]\n'.format(
                range_field,
                minval,
                maxval,
                range_unit))
    else:
        vrange_is_numeric = True

    # Sort vrange
    if vrange_is_numeric:
        vrange = tuple([
            min(minval, maxval),
            max(minval, maxval)])

    # Is vrange in practical range?
    vrange_min_is_practical = False
    if vrange_is_numeric:
        if vrange[0] < range_min:
            liner.send(
                '{}: {} to {} {} [MINIMUM VALUE SMALLER THAN {} {}]\n'.format(
                    range_field,
                    vrange[0],
                    vrange[1],
                    range_unit,
                    range_min,
                    range_unit))
        else:
            vrange_min_is_practical = True

    # Is maxtmelt in practical range?
    vrange_max_is_practical = False
    if vrange_min_is_practical:
        if vrange[1] > range_max:
            liner.send(
                '{}: {} to {} {} [MAXIMUM VALUE LARGER THAN {} {}]\n'.format(
                    range_field,
                    vrange[0],
                    vrange[1],
                    range_unit,
                    range_max,
                    range_unit))
        else:
            vrange_max_is_practical = True

    # Compute validity
    vrange_valid = all([
        vrange_is_numeric,
        vrange_min_is_practical,
        vrange_max_is_practical])

    # Show update
    if vrange_valid:
        liner.send(
            '{}: {} to {} {}\n'.format(
                range_field,
                vrange[0],
                vrange[1],
                range_unit))

    # Return validity
    return (
        vrange[0],
        vrange[1],
        vrange_valid)

def get_parsed_background(
    background,
    background_field,
    liner):
    '''
    Determine if background is valid.
    Internal use only.

    :: background
       type - string / db.vectorDB / None
       desc - path to background storage,
              or a vectorDB instance
    :: background_field
       type - string
       desc - background fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    background_field = _normalize_field(background_field)

    # Is background None?
    if background is None:
        liner.send(
            '{}: None Specified\n'.format(
                background_field))
        return True, None

    # Is background a vectorDB instance?
    if isinstance(background, db.vectorDB):
        liner.send(
            '{}: Contains {:,} Unique {}-mers\n'.format(
                background_field,
                len(background),
                background.K))
        return True, 'instance'

    # Is background a string?
    if isinstance(background, str):

        # Adjust background path
        indir = ut.get_adjusted_path(
            path=background,
            suffix='.oligopool.background')

        # Does background exist?
        background_exists = get_indir_validity(
            indir=indir,
            indir_suffix=None,
            indir_field=background_field,
            liner=liner)

        # Non-existent background
        if not background_exists:
            return False, None

        # Is background a valid vectorDB storage?
        else:

            # Open path as vectorDB instance
            try:
                vDB = db.vectorDB(
                    path=indir,
                    maximum_repeat_length=None)

            # Invalid attempt
            except Exception as E:
                liner.send(
                    '{}: {} [INVALID OR PRE-OPENED BACKGROUND OBJECT]\n'.format(
                        background_field,
                        _display_path(indir)))
                return False, None

            # Valid attempt
            else:
                liner.send(
                    '{}: Contains {:,} Unique {}-mers\n'.format(
                        background_field,
                        len(vDB),
                        vDB.K))
                vDB.close()
                return True, 'path'

    # Invalid input type
    else:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                background_field,
                background))
        return False, None

def get_parsed_backgrounds(
    backgrounds,
    backgrounds_field,
    liner):
    '''
    Validate single or multiple backgrounds.
    Returns (valid, info_list, max_K).
    Internal use only.

    :: backgrounds
       type - string / db.vectorDB / list / None
       desc - path(s) to background storage,
              or vectorDB instance(s), or list
              of paths/instances
    :: backgrounds_field
       type - string
       desc - backgrounds fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    backgrounds_field = _normalize_field(backgrounds_field)

    # Is backgrounds None?
    if backgrounds is None:
        liner.send(
            '{}: None Specified\n'.format(
                backgrounds_field))
        return True, None, None

    # Is backgrounds a single string or vectorDB?
    if isinstance(backgrounds, (str, db.vectorDB)):
        (valid, bg_type) = get_parsed_background(
            background=backgrounds,
            background_field=backgrounds_field,
            liner=liner)
        if not valid:
            return False, None, None
        # Get K value for single background
        if bg_type == 'instance':
            K = backgrounds.K
        elif bg_type == 'path':
            indir = ut.get_adjusted_path(
                path=backgrounds,
                suffix='.oligopool.background')
            vDB = db.vectorDB(path=indir, maximum_repeat_length=None)
            K = vDB.K
            vDB.close()
        else:
            K = None
        return True, [(backgrounds, bg_type, K)], K

    background_store = cx.deque()
    try:
        for bg in backgrounds:
            background_store.append(bg)
    except:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                backgrounds_field, backgrounds))
        return False, None, None

    # Empty list?
    if len(background_store) == 0:
        liner.send(
            '{}: None Specified\n'.format(
                backgrounds_field))
        return True, None, None

    # Single item in list? Delegate to single validation
    if len(background_store) == 1:
        return get_parsed_backgrounds(
            backgrounds=background_store[0],
            backgrounds_field=backgrounds_field,
            liner=liner)

    # Validate multiple backgrounds
    backgrounds_ok = True
    seen_paths = set()
    info_list = []
    max_K = 0

    # Prepare Spacing
    altspacing = ' '*len(backgrounds_field)

    # Show Header Update
    liner.send(
        '{}: {:,} Background Input(s)\n'.format(
            backgrounds_field, len(background_store)))

    # Core Validation Loop
    idx = 0
    while background_store:

        # Fetch background
        bg = background_store.popleft()
        idx += 1

        # Get path for duplicate detection
        if isinstance(bg, str):
            bg_path = ut.get_adjusted_path(
                path=bg,
                suffix='.oligopool.background')
        elif isinstance(bg, db.vectorDB):
            bg_path = bg.PATH.removesuffix(
                'vectorDB.ShareDB')
            bg_path = ut.removestarfix(
                string=bg_path,
                fix='/',
                loc=1)
        else:
            liner.send(
                '{}:  [{}] {} [INPUT TYPE IS INVALID]\n'.format(
                    altspacing, idx, bg))
            backgrounds_ok = False
            continue

        # Duplicate background?
        if bg_path in seen_paths:
            liner.send(
                '{}:  [{}] {} [DUPLICATE BACKGROUND]\n'.format(
                    altspacing, idx, _display_path(bg_path)))
            backgrounds_ok = False
            continue

        # Record path
        seen_paths.add(bg_path)

        # Validate this background
        if isinstance(bg, db.vectorDB):
            liner.send(
                '{}:  [{}] Contains {:,} Unique {}-mers\n'.format(
                    altspacing, idx, len(bg), bg.K))
            info_list.append((bg, 'instance', bg.K))
            max_K = max(max_K, bg.K)

        elif isinstance(bg, str):
            # Does background exist?
            background_exists = get_indir_validity(
                indir=bg_path,
                indir_suffix=None,
                indir_field='{}:  [{}]'.format(altspacing, idx),
                liner=liner)

            if not background_exists:
                backgrounds_ok = False
                continue

            # Is background a valid vectorDB storage?
            try:
                vDB = db.vectorDB(
                    path=bg_path,
                    maximum_repeat_length=None)
            except Exception as E:
                liner.send(
                    '{}:  [{}] {} [INVALID OR PRE-OPENED BACKGROUND OBJECT]\n'.format(
                        altspacing, idx, _display_path(bg_path)))
                backgrounds_ok = False
                continue
            else:
                liner.send(
                    '{}:  [{}] Contains {:,} Unique {}-mers\n'.format(
                        altspacing, idx, len(vDB), vDB.K))
                K = vDB.K
                vDB.close()
                info_list.append((bg, 'path', K))
                max_K = max(max_K, K)

    # Return Results
    if not backgrounds_ok:
        return False, None, None
    return True, info_list, max_K if max_K > 0 else None

def get_callback_validity(
    callback,
    callback_field,
    liner):
    '''
    Determine if callback function is valid.
    Internal use only.

    :: callback
       type - function / None
       desc - callback function specified
    :: background_field
       type - string
       desc - callback function fieldname
              used in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is callback None?
    if callback is None:
        liner.send(
            '{}: None Specified\n'.format(
                callback_field))
        return True

    # Is callback callable?
    if not callable(callback):
        liner.send(
            '{}: {} [INVALID FUNCTION]\n'.format(
                callback_field,
                callback))
        return False
    else:
        liner.send(
            '{}: Function w/ ID={}\n'.format(
                callback_field,
                id(callback)))
        return True

def get_errors_validity(
    errors,
    errors_field,
    errors_pre_desc,
    errors_post_desc,
    errors_base,
    indexfiles_valid,
    indexfiles,
    liner):
    '''
    Determine if numeric is a None or a
    positive Real. Internal use only.

    :: errors
       type - Real / None
       desc - error value to validate
    :: errors_field
       type - string
       desc - errors fieldname used in
              printing
    :: errors_pre_desc
       type - string
       desc - errors pre-description used
              in printing
    :: errors_post_desc
       type - string
       desc - errors post-description used
              in printing
    :: errors_base
       type - string
       desc - either 'A' or 'B'
    :: indexfile_valid
       type - boolean
       desc - True if indexfiles are valid,
              False otherwise
    :: indexfile
       type - iterable
       desc - filenames storing prepared
              indexes and models
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # errors is positive real?
    error_numeric = get_numeric_validity(
        numeric=errors,
        numeric_field=errors_field,
        numeric_pre_desc=errors_pre_desc,
        numeric_post_desc=errors_post_desc,
        minval=float('-inf'),
        maxval=float('+inf'),
        precheck=True,
        liner=liner)

    # errors status?
    if error_numeric:
        errors = round(errors)
        # errors is default
        if errors < 0.:
            if indexfiles_valid:
                errors = float('-inf')
                for indexfile in indexfiles:
                    archive = zf.ZipFile(
                        file=indexfile)
                    metamap = ut.loaddict(
                        archive=archive,
                        dfile='meta.map')
                    if errors_base == 'A':
                        errors = max(
                            errors,
                            metamap['associatetvalmax'])
                    else:
                        errors = max(
                            errors,
                            metamap['barcodetval'])
                    errors = round(errors)
                    archive.close()
                liner.send(
                    '{}:{}{}{} (Auto-Inferred)\n'.format(
                        errors_field,
                        errors_pre_desc,
                        errors,
                        errors_post_desc))
            else:
                liner.send(
                    '{}: Indeterminable\n'.format(
                        errors_field))
        # errors is specified
        else:
            liner.send('{}:{}{}{}\n'.format(
                errors_field,
                errors_pre_desc,
                errors,
                errors_post_desc))

    # Return errors validity
    return (errors,
        error_numeric)

def get_parsed_core_info(
    ncores,
    core_field,
    default,
    offset,
    liner):
    '''
    Determine if ncores is a positive Real and
    valid for given function. Internal use only.

    :: ncores
       type - Real
       desc - number of cores to use for function
    :: core_field
       type - string
       desc - ncores fieldname used in printing
    :: default
       type - None / integer
       desc - default value to use when ncores
              is Real but logically invalid
    :: offset
       type - integer
       desc - maximum numbers of core to withold
              from allocation if all cores usable
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # How many cores available in total?
    sys_cores = mp.cpu_count()

    # ncores is non-integer?
    if not isinstance(ncores, nu.Real):
        liner.send(
            '{}: Use {} out of {:,} Cores [INPUT TYPE IS INVALID]\n'.format(
                core_field, ncores, sys_cores))
        ncores_valid = False

    # ncores must be >= 0
    elif ncores < 0:
        liner.send(
            '{}: Use {:,} out of {:,} Cores [INPUT VALUE IS INVALID]\n'.format(
                core_field, ncores, sys_cores))
        ncores_valid = False

    # ncores is valid
    elif ncores >= 0:
        autoinferred = ''

        # ncores adjusted
        if  ncores > sys_cores:
            ncores = sys_cores
            autoinferred = '(Auto-Inferred)'

        # ncores defaulting
        elif ncores == 0:
            autoinferred = '(Auto-Inferred)'

            # Amdahl's Ironclad Law
            optncores = max(
                round(mt.sqrt(sys_cores)),
                round(mt.log2(sys_cores)))

            if default is None:
                ncores = optncores
            else:
                ncores = min(default, sys_cores)

        liner.send(
            '{}: Use {:,} out of {:,} Cores {}\n'.format(
                core_field,
                ncores,
                sys_cores,
                autoinferred))

        ncores = round(ncores)
        ncores_valid = True

    # Return adjusted ncores and validity
    adjcores = max(1, min(ncores, sys_cores-offset))
    return adjcores, ncores_valid

def get_parsed_memory_info(
    memlimit,
    memlimit_field,
    ncores,
    ncores_valid,
    liner):
    '''
    Determine if memlimit is a positive Real and
    valid for given function. Internal use only.

    :: memlimit
       type - Real
       desc - amount of memory to be allocated
              per core for a function
    :: core_field
       type - string
       desc - memlimit fieldname used in printing
    :: ncores
       type - integer
       desc - total number of cores used in the
              function
    :: ncores_valid
       type - boolean
       desc - if True ncores was parsed to be
              valid
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is ncores invalid?
    if not ncores_valid:
        liner.send(
            '{}: Use {} GB RAM per Core\n'.format(
                memlimit_field, memlimit))
        return memlimit, False

    # How much memory available in total?
    sys_mem = np.floor(
        pu.virtual_memory().total / (10**9)) - 2.0

    # memlimit is non-integer?
    if not isinstance(memlimit, nu.Real):
        liner.send(
            '{}: Use {} GB RAM per Core [INPUT TYPE IS INVALID]\n'.format(
                memlimit_field, memlimit))
        memlimit_valid = False

    # memlimit must be >= 0
    elif memlimit < 0:
        liner.send(
            '{}: Use {:.2f} GB RAM per Core [INPUT VALUE IS INVALID]\n'.format(
                memlimit_field, memlimit))
        memlimit_valid = False

    # memlimit is valid
    elif memlimit >= 0:
        autoinferred = ''

        # Normalize to Core Count
        sys_mem /= ncores

        # memlimit adjusted
        if  memlimit > sys_mem or \
            memlimit == 0.:
            memlimit = sys_mem
            autoinferred = '(Auto-Inferred)'

        liner.send(
            '{}: Use {:.2f} GB RAM per Core {}\n'.format(
                memlimit_field,
                memlimit,
                autoinferred))

        memlimit_valid = True

    # Return adjusted memlimit and validity
    return memlimit, memlimit_valid
