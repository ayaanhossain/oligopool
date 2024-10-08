'''
An example input for design_parser.py

Author: Ayaan Hossain
'''

import design_parser as dp
import pprint as pp

def get_promoter_list():
    '''Load core oligopool variants.'''
    with open('promoters.txt') as infile:
        promoter_list = [x.strip() for x in infile.readlines()]
    return promoter_list

def main():
    '''Design driver. Modify the spec here for design_parser execution.'''

    promoter_list = get_promoter_list()

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
                                       'type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 0,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer3',
                        'left_context_column': None,
                       'right_context_column': 'Cut1',
                            'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
            },

            'Cut1': {
                                     'type': 'motif',
                       'oligo_length_limit': 194,
                'motif_sequence_constraint': 'NNNGGATCCNNN',
                      'left_context_column': 'Primer1',
                     'right_context_column': 'Promoter',
                          'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
            },

            'Promoter': {
                     'type': 'variant',
                'sequences': promoter_list
            },

            'Barcode': {
                                    'type': 'barcode',
                      'oligo_length_limit': 250,
                          'barcode_length': 16,
                'minimum_hamming_distance': 3,
                   'maximum_repeat_length': 6,
                            'barcode_type': 1,
                         'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG'],
                     'left_context_column': 'Promoter',
                    'right_context_column': 'Primer2'
            },

            'Primer2': {
                                       'type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 0,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer3',
                        'left_context_column': 'Barcode',
                       'right_context_column': 'Cut2',
                            'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
            },

            'Cut2': {
                                     'type': 'motif',
                       'oligo_length_limit': 194,
                'motif_sequence_constraint': 'NNNTCTAGANNN',
                      'left_context_column': 'Primer2',
                     'right_context_column': 'Primer3',
                          'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
            },

            'Primer3': {
                                       'type': 'primer',
                         'oligo_length_limit': 250,
                 'primer_sequence_constraint': 'NNNNNNNNNNNNNNNNNNNNNN',
                                'primer_type': 1,
                'minimum_melting_temperature': 53,
                'maximum_melting_temperature': 55,
                      'maximum_repeat_length': 10,
                       'paired_primer_column': 'Primer2',
                        'left_context_column': 'Cut2',
                       'right_context_column': 'Filler',
                            'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
            },

            'Filler': {
                                'type': 'spacer',
                  'oligo_length_limit': 250,
                       'spacer_length': None,
                 'left_context_column': 'Primer3',
                'right_context_column': None,
                     'excluded_motifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']
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

    # for df in output:
    #     print(output[df])

if __name__ == '__main__':
    main()