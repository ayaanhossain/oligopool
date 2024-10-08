'''
An Example of Oligopool Calculator Automated Design Mode Parser.

Use the `design_parser(...)` method for pipeline execution.

Look into run_design_parser.py for an example of
how to use this parser.

Author: Ayaan Hossain
'''

import os
import shutil
import uuid
import atexit
import contextlib

import collections
import pandas

import oligopool as op


# --= WORKSPACE METHODS =--

@contextlib.contextmanager
def ignored(*exceptions):
    '''Raymond Hettinger's ignored context.'''
    try:
        yield
    except exceptions:
        pass

def removestarfix(string, fix, loc):
    '''Remove prefix or suffix from string.'''
    if loc <= 0:
        if string.startswith(fix):
            return string[len(fix):]
    elif loc > 0:
        if string.endswith(fix):
            return string[:-len(fix)]
    return string

def get_adjusted_path(path, suffix):
    '''Adjust output_file with suffix.'''
    if not isinstance(path, str):
        return None
    path = path.strip()
    path = removestarfix(
        string=path,
        fix='/',
        loc=1)
    if not suffix is None and \
    not path.endswith(suffix):
        path += str(suffix)
    return path

def remove_directory(dirpath):
    '''Delete dirpath and its contents.'''
    if not dirpath is None:
        dirpath = get_adjusted_path(
            path=dirpath,
            suffix=None)
        with ignored(OSError):
            shutil.rmtree(dirpath)

def setup_directory(dirpath):
    '''Build dirpath.'''
    dirpath = get_adjusted_path(
        path=dirpath,
        suffix=None)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def setup_workspace(workspace_name):
    '''Setup temporary workspace.'''
    atexit.register(
        remove_directory,
        workspace_name)
    setup_directory(
        dirpath=workspace_name)


# --= DESIGN MODE METHODS =--

def setup_dataframe(pool_size, element_names):
    '''Build the initial DataFrame with placeholders.'''

    # Blank DataFrame
    dataframe = pandas.DataFrame()

    # Build IDs
    dataframe.insert(
        0, 'ID',  ['V-{}'.format(idx+1) for idx in range(pool_size)])

    # Build Element Columns
    for idx,colname in enumerate(element_names):
        dataframe.insert(idx, colname, '-')

    # Re-Index
    dataframe.set_index('ID', inplace=True)

    # Return Results
    return dataframe

def setup_background(workspace_name, background_spec, stats_dict):
    '''Build the background within workspace.'''

    # Build Background
    status = False
    if (not (background_spec is None)) and \
       (len(background_spec) > 0):

        # Define Background Name
        background_name = workspace_name + '/oligopool_background'

        # Build Background
        try:
            stats = op.background(
                input_data=background_spec['input_data'],
                maximum_repeat_length=background_spec['maximum_repeat_length'],
                output_directory=background_name)
            stats_dict['background'] = stats
            status = True
        except:
            stats = {
                'status'   : False,
                'basis'    : 'infeasible',
                'step'     : 0,
                'step_name': 'parsing-input',
                'vars'     : {},
                'warns'    : {}}
            status = False
    else:
        print()
        background_name = None

    # Return Results
    return status, background_name, stats_dict

def insert_variants(dataframe, elements_spec, stats_dict):
    '''Insert variants into DataFrame.'''

    # Process All Variants
    status = False
    for element_name in elements_spec:
        if elements_spec[element_name]['element_type'] == 'variant':
            try:
                dataframe[element_name] = elements_spec[element_name]['sequences']
                stats_dict['oligopool'][element_name] = {
                    'status'   : True,
                    'basis'    : 'solved',
                    'step'     : 1,
                    'step_name': 'inserting-variants',
                    'vars'     : {},
                    'warns'    : {}
                }
                status = True
            except:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'   : False,
                    'basis'    : 'infeasible',
                    'step'     : 0,
                    'step_name': 'parsing-input',
                    'vars'     : {},
                    'warns'    : {}
                }
                status = False
                break

    # Show Update
    print('[Oligopool Calculator: Design Mode - Inserting Variants]\n')
    print(dataframe)
    print()

    # Return Results
    return status, dataframe, stats_dict

def insert_motifs(dataframe, elements_spec, stats_dict):
    '''Insert any degenerate motifs into DataFrame.'''

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['element_type'] == 'motif':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = op.motif(
                    input_data=dataframe,
                    oligo_length_limit=elements_spec[element_name]['oligo_length_limit'],
                    motif_sequence_constraint=elements_spec[element_name]['motif_sequence_constraint'],
                    motif_column=element_name,
                    output_file=None,
                    left_context_column=elements_spec[element_name]['left_context_column'],
                    right_context_column=elements_spec[element_name]['right_context_column'],
                    excluded_motifs=elements_spec[element_name]['excluded_motifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except Exception as E:
                # raise E
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'   : False,
                    'basis'    : 'infeasible',
                    'step'     : 0,
                    'step_name': 'parsing-input',
                    'vars'     : {},
                    'warns'    : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def get_primer_order(elements_spec):
    '''What is the optimal order of designing the primers?'''
    order = collections.deque()
    for element_name in elements_spec:
        if elements_spec[element_name]['element_type'] == 'primer':
            if not element_name in order:
                order.append(element_name)
            paired_primer = elements_spec[element_name]['paired_primer_column']
            if not paired_primer is None and \
               not paired_primer in order:
                index = order.index(element_name)
                order.insert(index, paired_primer)
    return order

def insert_primers(element_names, background, dataframe, elements_spec, stats_dict):
    '''Insert primers into DataFrame.'''

    # Process All Variants
    status = False
    incomplete = False
    order = get_primer_order(elements_spec=elements_spec)
    # print(order)
    # print(dataframe.columns)
    # print(element_names)
    while order:
        element_name = order.popleft()
        try:
            if element_name in dataframe.columns:
                dataframe.drop(element_name, inplace=True, axis=1)
            paired_primer_column = elements_spec[element_name]['paired_primer_column']
            if paired_primer_column in order:
                paired_primer_column = None # element_name is a Predecessor
            # Left Context is Sensitive
            if elements_spec[element_name]['left_context_column'] in stats_dict:
                left_context_column = elements_spec[element_name]['left_context_column']
            else:
                left_context_column = None
            # Right Context is Sensitive
            if elements_spec[element_name]['right_context_column'] in stats_dict:
                right_context_column = elements_spec[element_name]['right_context_column']
            else:
                right_context_column = None
            dataframe, stats = op.primer(
                input_data=dataframe,
                oligo_length_limit=elements_spec[element_name]['oligo_length_limit'],
                primer_sequence_constraint=elements_spec[element_name]['primer_sequence_constraint'],
                primer_type=elements_spec[element_name]['primer_type'],
                minimum_melting_temperature=elements_spec[element_name]['minimum_melting_temperature'],
                maximum_melting_temperature=elements_spec[element_name]['maximum_melting_temperature'],
                maximum_repeat_length=elements_spec[element_name]['maximum_repeat_length'],
                primer_column=element_name,
                output_file=None,
                paired_primer_column=paired_primer_column,
                left_context_column=left_context_column,
                right_context_column=right_context_column,
                excluded_motifs=elements_spec[element_name]['excluded_motifs'],
                background_directory=background,
                verbose=True)
            # print(dataframe)
            # print(dataframe.columns)
            for col in element_names:
                if not col in dataframe.columns:
                    # print(col)
                    dataframe[col] = '-'
            # print(dataframe.columns)
            dataframe = dataframe[element_names]
            # print(dataframe)
            stats_dict['oligopool'][element_name] = stats
            status = True
            if stats['status'] is False:
                incomplete = True
                break
        except Exception as E:
            # raise E
            dataframe = None
            stats_dict['oligopool'][element_name] = {
                'status'   : False,
                'basis'    : 'infeasible',
                'step'     : 0,
                'step_name': 'parsing-input',
                'vars'     : {},
                'warns'    : {}}
            status = False
            break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def insert_barcodes(dataframe, elements_spec, stats_dict):
    '''Insert barcodes into DataFrame.'''

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['element_type'] == 'barcode':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = op.barcode(
                    input_data=dataframe,
                    oligo_length_limit=elements_spec[element_name]['oligo_length_limit'],
                    barcode_length=elements_spec[element_name]['barcode_length'],
                    minimum_hamming_distance=elements_spec[element_name]['minimum_hamming_distance'],
                    maximum_repeat_length=elements_spec[element_name]['maximum_repeat_length'],
                    barcode_column=element_name,
                    output_file=None,
                    barcode_type=elements_spec[element_name]['barcode_type'],
                    left_context_column=elements_spec[element_name]['left_context_column'],
                    right_context_column=elements_spec[element_name]['right_context_column'],
                    excluded_motifs=elements_spec[element_name]['excluded_motifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except Exception as E:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 0,
                    'step_name': 'parsing-input',
                    'vars'    : {},
                    'warns'   : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def insert_spacers(dataframe, elements_spec, stats_dict):
    '''Insert spacers into DataFrame.'''

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['element_type'] == 'spacer':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = op.spacer(
                    input_data=dataframe,
                    oligo_length_limit=elements_spec[element_name]['oligo_length_limit'],
                    spacer_column=element_name,
                    output_file=None,
                    spacer_length=elements_spec[element_name]['spacer_length'],
                    left_context_column=elements_spec[element_name]['left_context_column'],
                    right_context_column=elements_spec[element_name]['right_context_column'],
                    excluded_motifs=elements_spec[element_name]['excluded_motifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'   : False,
                    'basis'    : 'infeasible',
                    'step'     : 0,
                    'step_name': 'parsing-input',
                    'vars'     : {},
                    'warns'    : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def split_oligos(dataframe, split_spec, stats_dict):
    '''Split long oligos into shorter fragments.'''

    # Process All Variants
    status = False
    incomplete = False
    try:
        dataframe, stats = op.split(
            input_data=dataframe,
            split_length_limit=split_spec['split_length_limit'],
            minimum_melting_temperature=split_spec['minimum_melting_temperature'],
            minimum_hamming_distance=split_spec['minimum_hamming_distance'],
            minimum_overlap_length=split_spec['minimum_overlap_length'],
            maximum_overlap_length=split_spec['maximum_overlap_length'],
            output_file=None)
        stats_dict['split'] = stats
        status = True
        if stats['status'] is False:
            incomplete = True
    except:
        dataframe = None
        stats_dict['split'] = {
            'status'   : False,
            'basis'    : 'infeasible',
            'step'     : 0,
            'step_name': 'parsing-input',
            'vars'     : {},
            'warns'    : {}}
        status = False

    # Return Results
    return status, incomplete, dataframe, stats_dict

def pad_oligos(dataframe, padding_spec, stats_dict):
    '''Add padding to split oligos.'''

    # Process All Variants
    status = False
    incomplete = False
    paddedframes = []
    stats_dict['padding'] = {}
    for split_column in dataframe.columns:
        try:
            paddedframe, stats = op.pad(
                input_data=dataframe,
                split_column=split_column,
                typeIIS_system=padding_spec['typeIIS_system'],
                oligo_length_limit=padding_spec['oligo_length_limit'],
                minimum_melting_temperature=padding_spec['minimum_melting_temperature'],
                maximum_melting_temperature=padding_spec['maximum_melting_temperature'],
                maximum_repeat_length=padding_spec['maximum_repeat_length'],
                output_file=None)
            stats_dict['padding'][split_column] = stats
            status = True
            if stats['status'] is False:
                incomplete = True
                break
        except:
            paddedframe = None
            stats_dict['padding'][split_column] = {
                'status'   : False,
                'basis'    : 'infeasible',
                'step'     : 0,
                'step_name': 'parsing-input',
                'vars'     : {},
                'warns'    : {}}
            status = False
            break
        paddedframes.append(paddedframe)

    # Return Results
    return status, incomplete, paddedframes, stats_dict

def design_parser(
    pool_size:int,
    element_names:list,
    elements_spec:dict,
    background_spec:dict,
    split_spec:dict,
    padding_spec:dict) -> dict:
    '''Oligopool Calculator Design Mode Automated Pipeline.

    This function parses design specifications for an oligopool and processes them
    to generate the desired library.

    Parameters:
    -----------
        pool_size : int
            The number of unique variants in the oligopool.
            Example: pool_size = 4350

        element_names : list
            A list of names for each oligopool element in the order of desired adjacency.
            Example: element_names = ['promoter', 'gene', 'terminator']

        elements_spec : dict
            A dictionary containing various element design parameters for each named element.
            The keys should correspond to the names in element_names.
            Example: elements_spec = {'primer': {...}, 'promoter': {...}, 'barcode': {...}, ...}

        background_spec : dict
            A dictionary with background parameters for the oligopool design.

        split_spec : dict
            A dictionary specifying how to split or segment the oligopool design.

        padding_spec : dict
            A dictionary containing specifications for padding sequences in the design.


    Returns:
    --------
        A dictionary containing all results and pipeline statistics in a dictionary.


    Notes:
    ------
        See 'run_design_parser.py' for example for an input structure.
    '''

    # Book-keeping
    stats_dict = {
        'background': None,
         'oligopool': {},
             'split': None,
           'padding': None}
    output = None

    # Setup Temporary Workspace
    workspace_name = str(uuid.uuid4())
    setup_workspace(
        workspace_name=workspace_name)

    # Build Initial DataFrame
    dataframe = setup_dataframe(
        pool_size=pool_size,
        element_names=element_names)
    print('\n[Oligopool Calculator: Design Mode - Initialized DataFrame]\n')
    print(dataframe)
    print()

    # Setup Background
    if ((background_spec is None) or \
        (len(background_spec) <= 0)):
        background_name = None
        print()
    else:
        (status,
        background_name,
        stats_dict) = setup_background(
            workspace_name=workspace_name,
            background_spec=background_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                'step': 1,
                'desc': 'Step 1: Background Input Error',
                'output': output,
            'stats_dict': stats_dict}
        elif stats_dict['background']['status'] is False:
            return {
                'step': 1,
            'step_name': 'Step 1: Background Unsuccessful',
                'output': output,
            'stats_dict': stats_dict}

    # Insert Core Variants
    (status,
    dataframe,
    stats_dict) = insert_variants(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 2,
         'step_name': 'Step 2: Variants Input Error',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Motifs
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_motifs(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 3,
         'step_name': 'Step 3: Motif Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 3,
         'step_name': 'Step 3: Motif Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Primers
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_primers(
        element_names=element_names,
        background=background_name,
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 4,
         'step_name': 'Step 4: Primer Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 4,
         'step_name': 'Step 4: Primer Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Barcodes
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_barcodes(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 5,
         'step_name': 'Step 5: Barcode Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 5,
         'step_name': 'Step 5: Barcode Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Spacers
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_spacers(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 6,
         'step_name': 'Step 6: Spacer Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 6,
         'step_name': 'Step 6: Spacer Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Finalize DataFrame
    final_dataframe, _ = op.final(
        input_data=dataframe,
        output_file=None,
        verbose=True)
    final_dataframe.drop('OligoLength', inplace=True, axis=1)

    # Split Oligos
    split_dataframe = None
    if (not (split_spec is None)) or \
       (len(split_spec) > 0):
        (status,
        incomplete,
        split_dataframe,
        stats_dict) = split_oligos(
            dataframe=final_dataframe,
            split_spec=split_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                  'step': 8,
             'step_name': 'Step 8: Split Input Error',
                'output': output,
            'stats_dict': stats_dict}
        if incomplete:
            return {
                  'step': 8,
             'step_name': 'Step 8: Split Infeasible',
                'output': output,
            'stats_dict': stats_dict}

    # Pad Oligos
    padded_dataframes = None
    if ((not (split_spec is None)) and \
        (len(split_spec) > 0)) and \
       ((not (padding_spec is None)) and \
        (len(padding_spec) > 0)):
        (status,
        incomplete,
        padded_dataframes,
        stats_dict) = pad_oligos(
            dataframe=split_dataframe,
            padding_spec=padding_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                  'step': 9,
             'step_name': 'Step 9: Padding Input Error',
                'output': output,
            'stats_dict': stats_dict}
        if incomplete:
            return {
                  'step': 9,
             'step_name': 'Step 9: Padding Infeasible',
                'output': output,
            'stats_dict': stats_dict}

    # Remove Workspace
    remove_directory(dirpath=workspace_name)

    # Compute Ouput
    output = {
        'annotated_complete': dataframe.to_dict(
            orient='index'),
         'presplit_complete': final_dataframe.to_dict(
             orient='index'),
           'padded_complete': [
               df.to_dict(orient='index') for df in padded_dataframes]
    }

    # Return Results
    return {
              'step': 10,
         'step_name': 'Step 10: Design Infeasible',
            'output': output,
        'stats_dict': stats_dict}


# --= MISC TEST =--

def test_primer_order():
    '''Resolve the dependency graph of a system of primers.'''

    # Define all Primers
    elements_spec = {
        'Primer1': {
            'element_type': 'primer',
            'paired_primer_column': 'Primer2',
        },

        'Primer2': {
            'element_type': 'primer',
            'paired_primer_column': 'Primer1',
        },

        'Primer3': {
            'element_type': 'primer',
            'paired_primer_column': 'Primer2',
        },
    }

    # Resolve and show order
    order = get_primer_order(elements_spec)
    print(order)

if __name__ == '__main__':
    # test_primer_order()
    help(design_parser)