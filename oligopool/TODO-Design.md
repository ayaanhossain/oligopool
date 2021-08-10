codemeta
~~* edge effect refactor to utils.py?~~
~~* Checking constraints in Maker switchable~~
~~* Extract test for modules in tests directory~~
  ~~* barcode.py~~
  ~~* primer.py~~
~~* Failure states are recorded in dictionary
  - type: infeasible, unsolved, solved
  - step: 0  for solved,
          -1 for unsolved,
          xx for infeasible step integer id
  - variables: {
      keys: values
    } // empty when successful~~
~~* Column Requirement specification for parsing~~
~~* introduce search order, an integer paramter
  which is 10^exponent of failure count / retry attempts~~
~~* exmotifs come after context~~
~~* Add everything before and after context as part of elongated context
  -- e.g. if LeftContext is ColA, then everything before and including
     ColA is left context; this builds full left context~~
~~* Introduce '-' as a valid DNA base signifying a gap base, which
  is used to denote zero length spacers, and removed once oligopool
  finalized~~

motif      -- done
background -- done
primer     -- done
barcode    -- done
spacer     -- done
split      -- done
pad        -- done
lenstat    -- done
final      -- done

Refactoring
~~* Build exmotif partition and inverse partition into
  parsing of exmotifs itself, based on left and right
  context inputs (treated as binary presence)~~
~~* Parsing exmotifs then take the partitions as input
  instead of building them in parsing edgeeffects~~
~~* introduce a merge boolean flag to merge all prefix
  and suffix dict entries into prefix and suffixsets
  which will be used in primer design objectives~~

WARNINGS
~~* Introduce stepname with every stats dictionary~~
~~* Explicitly add step number to pipeline sections~~
~~* if any left paritions are right blocked or any right
  partitions are left blocked then warn from parsing
  of exmotifs~~
~~* show impact of edgeeffects count (already done), but
  also report the 5' and 3' constants from sequence
  constraints that create those edgeeffects~~

motif.py
~~add a fixed length motif or constant sequence
between two columns or at the 3' end by default
default 3' end~~

spacer.py
~~add variable length spacer between two cols of
choice or at the ends of the oligos, preventing
edge effects and upto a fixed length oligo~~

finalize.py
~~concatenate all cols, and run length checks
for oligopool length bounds~~

barcode.py
~~* Revise barcode.py~~
~~* Find a circular queue alternative to bitarray approach for assignment~~
~~* Return all failures broken down from engine~~
~~* Single-point Jumper Computation~~
~~* Check Barcode Design Feasibility~~
  ~~- if barcode length is infeasible for targetsize~~
~~* If both left and right contexts are constant then break after 1 check only,~~
  ~~for other cases just break after all carr checks~~
~~* Check if context has exmotifs in them?~~
  ~~-- Better approach, take max context len~~
     ~~and then approve barcode if exmotif~~
     ~~before or beyond barcode (doesn't touch~~
     ~~the barcode either in starting or ending)~~
~~* Refactor all related functions~~
~~* Build out final module interface~~
  ~~-- doc string~~
  ~~-- value parsing~~
  ~~-- engine processing~~
  ~~-- output dumping / return~~
~~* Build out statistics in barcode.py~~
~~* Report Time Elapsed at Construction End~~

primer.py
~~* Start developing primer.py~~
~~* Determine if pool-universal primer possible for edge effects~~
~~* Determine if Tm feassible primer possible given sequence constraint~~
~~* Finish motif-edge effect (context) objective function~~
~~* Finish Tm objective function~~
~~* Finish motif objective function~~
~~* Update NRPCalc Fold to include Cofold~~
~~* Finish homodimer objective function~~
~~* Finish co-fold objective function~~
~~* Check if sequence constraint has exmotifs in them?~~
~~* Finish interface~~
~~* Parse Exmotifs, then Context, then Edge-Effect (requires both fulfilled)~~
~~* Eliminate Primer Lmax from parameter (use 4 or 5)?~~
  ~~-- background has its own Lmax~~
  ~~-- auto-infer based on longest common ATGC string between~~
     onstraints?~~
     ~~-- remove all characters non-ATGC~~
     ~~-- align every substring against everything
        else and it's reverse complement~~
     ~~-- longest common substring will be found~~
     ~~-- Lmax = this length~~
~~* RNAduplex to check for energy of the longest contig (Alex strat)~~
  ~~-- introduce functionality in NRP Calculator~~
  ~~-- >> NO NEED TO DO THIS~~
* Move over Alex's changes to shared PyVRNA to PyVRNA directory
~~* Check if sequence constraint is degenerate enough?~~
  ~~-- maybe this and other type checking becomes part of valparse.py functions?~~
~~* Ensure all fail states and causes recorded~~
~~* Polish out primer.py~~
~~* Check if context has exmotifs in them?~~
~~* Build out statistics~~
~~* Interface has left and right context names, and the option
  for a FASTA file, or a single DNA string for vector~~
~~* Primers non-repetitive to everything before and after Left and Right Context~~
~~* Report Time Elapsed at Construction End
~~* ignore palindromic hexamers, and alter optimization criteria
  if cutsites are part of primer sequence constraint~~

primer - padding common
~~* Refactor all core primer functions into coreprime.py~~
~~* Return not just the primer/pad sequence, but the insert~~

padding.py
~~* Start developing padding.py~~
~~* Build out type2s dictionary~~
~~* Pad constraints must have 6 Ns on the left or right depending on type~~
~~* Build out the Local Model Function~~
~~* Build and test Forward Pad~~
~~* Build and test Reverse Pad~~
~~* Follow Gameplan~~

split.py
~~* Pad shorter sequences, throwing warnings~~
    ~~- perhaps a validation?~~
    ~~- perhaps this is dynamic and determined
      as each split is made; a split comprising
      of less than T nucleotides of actual
      fragment content is flagged out~~
~~* Chop off padding before returning~~