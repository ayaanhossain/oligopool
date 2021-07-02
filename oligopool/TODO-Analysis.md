codemeta
* All stats are output as dictionary
* Optionally engines may compute stats dict

pack.py

* Two Types
    - Store -- Store R1 and revcomp(R2) as a tuple
    - Merge -- Merge R1 and revcomp(R2) via SHW alignment
               then store the merged read as a tuple

* Merging is two types:
    Both attempted in Union > Intersection order
    Union and Intersection Merge
    Ref: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2579-2/figures/1

* Option to subselect reads based on LeftConstant and RightConstant

index.py (same as prep.py)

~~* Indexing Structure
    - Each index takes in barcode.csv and variant.csv
    - ID must be shared between barcode and variants
    - Fields: ID LeftConstant Barcode(/Variant) RightConstant
    - Either LeftConstant or RightConstant must be present~~

~~* Introduce get_paired_parsed_indata_info(...) as a superset of
  get_parsed_indata_info(...) that ensures paired df indexes
  match each other and reindex the input df based on paired
  df. Look into get_parsed_spacerlen_info(...) for clues.~~

~~* Finish Engine~~
~~* Finish Scry~~

count.py

* Counting is of three types:
    - Union Counting: Map and count all variants for each read,
      applied separately -- reads contain elements of one or more
      indices applied
    - Intersection Counting: Map and count all variants jointly
      for each read, applied jointly (tuple counting) -- reads
      contain elements of all indices applied
    - Exclusive Counting: Map and count variants exclusively for
      each read, applied sequentially assigning read to one index
      whichever matches best -- reads expected to contain elements
      from exactly one index (e.g. Spike-In or Mixed Sequencing)

* Unlike Dx-Seq, there is no input of errors for the barcodes to be
  counted -- it is inferred during indexing, and becomes fixed

* Outputs
    - Union - one CSV file for each index
    - Intersection - one CSV with with a tuple of indices

