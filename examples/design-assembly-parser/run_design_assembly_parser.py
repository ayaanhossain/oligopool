'''
An example input for design_assembly_parser.py

Author: Ayaan Hossain
'''

import design_assembly_parser as dp
import pprint as pp

def get_promoter_list():
    '''Load core promoter variants.'''
    with open('promoters.txt') as infile:
        promoter_list = [x.strip() for x in infile.readlines()]
    return promoter_list

def run_design_parser():
    '''Design driver. Modify the spec here for design_parser execution.'''

    promoter_list = get_promoter_list()   # Load Core Oligo Variants

    excluded_motifs = ['GGATCC', 'TCTAGA', 'GGTCTC', 'GAGACC', 'CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
    promoter_list = [
        p for p in promoter_list
        if not any(e in p for e in excluded_motifs)
    ]

    # REMEMBER TO PROCESS OR WRITE THE OUTPUT from the line below

    output = dp.design_parser(

        pool_size=len(promoter_list),

        element_names=[
            'Primer1',
            'Cut1',
            'Promoter',
            'Barcode',
            'Primer2',
            'Cut2',
            'Primer3',
            'Filler'],

        elements_spec={
            'Primer1': {
                               'element_type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 0,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer3',
                        'left_context_column': None,
                       'right_context_column': 'Cut1',
                            'excluded_motifs': excluded_motifs
            },

            'Cut1': {
                             'element_type': 'motif',
                       'oligo_length_limit': 194,
                'motif_sequence_constraint': 'NNNGGATCCNNN',
                    'maximum_repeat_length': 8,
                      'left_context_column': 'Primer1',
                     'right_context_column': 'Promoter',
                          'excluded_motifs': excluded_motifs
            },

            'Promoter': {
                     'element_type': 'variant',
                'sequences': promoter_list         # Your Core Oligo Variants go here
            },

            'Barcode': {
                            'element_type': 'barcode',
                      'oligo_length_limit': 250,
                          'barcode_length': 16,
                'minimum_hamming_distance': 3,
                   'maximum_repeat_length': 8,
                            'barcode_type': 1,
                         'excluded_motifs': excluded_motifs,
                     'left_context_column': 'Promoter',
                    'right_context_column': 'Primer2'
            },

            'Primer2': {
                               'element_type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 0,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer3',
                        'left_context_column': 'Barcode',
                       'right_context_column': 'Cut2',
                            'excluded_motifs': excluded_motifs
            },

            'Cut2': {
                             'element_type': 'motif',
                       'oligo_length_limit': 194,
                'motif_sequence_constraint': 'NNNTCTAGANNN',
                    'maximum_repeat_length': 8,
                      'left_context_column': 'Primer2',
                     'right_context_column': 'Primer3',
                          'excluded_motifs': excluded_motifs
            },

            'Primer3': {
                               'element_type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 1,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer2',
                        'left_context_column': 'Cut2',
                       'right_context_column': 'Filler',
                            'excluded_motifs': excluded_motifs
            },

            'Filler': {
                         'element_type': 'spacer',
                   'oligo_length_limit': 250,
                'maximum_repeat_length': 8,
                        'spacer_length': None,
                  'left_context_column': 'Primer3',
                 'right_context_column': None,
                      'excluded_motifs': excluded_motifs
            },
        },

        background_spec={
                       'input_data': promoter_list,
            'maximum_repeat_length': 12
        },

        split_spec={
                     'split_length_limit': 170,
            'minimum_melting_temperature': 40,
               'minimum_hamming_distance': 1,
                 'minimum_overlap_length': 20,
                 'maximum_overlap_length': 50
        },

        padding_spec={
                         'typeIIS_system': 'bsaI',
                     'oligo_length_limit': 200,
            'minimum_melting_temperature': 40,
            'maximum_melting_temperature': 80,
                  'maximum_repeat_length': 15,
        }
    )

    pp.pprint(output['stats_dict'])

    # REMEMBER TO PROCESS OR WRITE THE OUTPUT from design_parser

if __name__ == '__main__':
    run_design_parser()