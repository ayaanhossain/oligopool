import os

import oligopool as op


def main():
    '''
    Example analysis pipeline showing Combinatorial Counting.

    For details, use help(...) over the methods or check out
    our read_simulator.ipynb and OligopoolCalculatorInAction.ipynb
    notebooks describing the data generation and analysis logic.
    '''

    '''
    Step 1: Index all Barcodes to be Mapped
    '''

    # Index the First Barcode
    op.index(
        barcode_data='ribozyme_architecture.csv',    # Filename storing barcode information
        barcode_column='BC1',                        # Index BC1 only
        barcode_prefix_column='OrangeForwardPrimer', # which should be next to OrangeForwardPrimer
        barcode_suffix_column=None,                  # We have no flanking right constant for BC1
        barcode_prefix_gap=0,                        # BC1 should be right next to prefix
        barcode_suffix_gap=0,
        index_file='BC1',                            # Store results in BC1.oligopool.index file
        verbose=True,
    )

    # Index the Second Barcode
    op.index(
        barcode_data='ribozyme_architecture.csv',
        barcode_column='BC2',                        # Index BC2 only
        barcode_prefix_column='PinkForwardPrimer',   # which should be next to PinkForwardPrimer
        barcode_suffix_column='YellowReversePrimer', # and YellowReversePrimer as flanking right constant
        barcode_prefix_gap=0,
        barcode_suffix_gap=0,
        associate_data='ribozyme_architecture.csv',  # Associates are in same DataFrame
        associate_column='Variant',
        associate_prefix_column=None,
        associate_suffix_column='PinkForwardPrimer',
        index_file='BC2',                            # Store results in BC2.oligopool.index file
    )

    '''
    Step 2: Pack all FASTQ Files to be Counted
    '''

    # If we have multiple pairs of FASTQ files, execute this in a loop
    # for each pair serially
    op.pack(
        r1_fastq_file='ribozyme_1M_R1.fq.gz',
        r2_fastq_file='ribozyme_1M_R2.fq.gz',
        r1_read_type=0, # R1 is in Forward Orientation
        r2_read_type=1, # R2 is in Reverse Orientation
        minimum_r1_read_quality=30, # Filter out any reads with Phred Score less than 30
        minimum_r2_read_quality=30,
        minimum_r1_read_length=10,  # Filter out any reads shorter than 10 bases
        minimum_r2_read_length=10,
        pack_type=1,   # R1 and R2 needs to be merged
        pack_size=0.1, # Store 100K reads per read pack
        pack_file='NGS',
        core_count=4,  # Use only 4 CPU cores
    )

    '''
    Step 3: Count all Packed Reads (Step 2) using Index (Step 1)
    '''

    # Custom Analysis Function
    def myfunc(read, ID, count, coreid):
        '''
        A custom read processing function that is
        executed concurrently with counting, which
        performs additional analysis, and returns
        a boolean indicating whether a read should
        be accepted or not.

        These functions are useful when one needs
        to extract additional information from the
        reads being counted, or specify additional
        criteria for accepting a read for counting.

        All custom counting functions must at least
        accept the following arguments.

        :: read
           type - string
           desc - the read being counted / analyzed
        :: ID
           type - tuple
           desc - a tuple of IDs, one for each element
                  mapped from each of the given indexes
                  e.g. ('Index-1-Barcode-47',
                        'Index-2-Barcode-12',
                        '-',
                        'Index-4-Barcode-77)
                  is a potential ID for a given read
                  when four indexes are specified for
                  counting, and none of the elements
                  from the third index got mapped.
                  The '-' indicates a missing value.
        :: count
           type - integer
           desc - the associated read / ID frequency
                  count in a given read pack
        :: coreid
           type - integer
           desc - the coreid integer identifier that
                  processed this read
        '''
        return True # Accept everything, for now.

    # The Count Matrix will be saved to 'CC.oligopool.xcount.csv'
    op.xcount(
        index_files=['BC1', 'BC2'],     # Map combinations of BC1 and BC2
        pack_file='NGS.oligopool.pack',
        count_file='CC',
        mapping_type=1,
        barcode_errors=-1,
        callback=None,
        core_count=0,
        memory_limit=0.0,
        verbose=True)

    # Cleanups
    os.remove('BC1.oligopool.index')
    os.remove('BC2.oligopool.index')
    os.remove('NGS.oligopool.pack')
    os.remove('CC.oligopool.xcount.csv')


if __name__ == '__main__':
    main()
