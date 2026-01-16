#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import json
import argparse
import textwrap
import functools
import datetime as dt

from . import __version__, __author__
from .background import background
from .barcode import barcode
from .primer import primer
from .motif import motif
from .spacer import spacer
from .merge import merge
from .revcomp import revcomp
from .lenstat import lenstat
from .final import final
from .split import split
from .pad import pad
from .index import index
from .pack import pack
from .acount import acount
from .xcount import xcount


__doc__ = f'oligopool v{__version__}\nby ah'

CLI_MANUAL = '''
Oligopool CLI Manual

Usage:
  oligopool COMMAND --argument=<value> ...
  op COMMAND --argument=<value> ...

Notes:
  - Design/transform commands require --output-file in CLI mode.
  - Run "oligopool COMMAND" to see command-specific options.
  - Use "oligopool manual topics" to list manual topics.
'''

# Banner toggles to avoid duplicate prints in embedded contexts.
header = False
footer = False


class OligopoolParser(argparse.ArgumentParser):
    '''
    Custom Parser to show help messages on
    error. Internal use only.
    '''
    def error(self, message):
        if message:
            message = message.strip()
            if not message.endswith(':'):
                print(f'error:\n{message}\n')
        try:
            self.print_help()
            print()
        except Exception:
            pass
        sys.exit(404)

    def parse_known_args(self, args=None, namespace=None):
        """Treat unknown args as a hard error for subparsers."""
        args, argv = super().parse_known_args(args, namespace)
        if argv:
            argv_string = ' '.join(argv)
            if argv_string:
                self.error(f'unrecognized arguments: {argv_string}')
        return args, argv


class OligopoolFormatter(argparse.RawTextHelpFormatter):
    '''
    Custom HelpFormatter for OligopoolParser.
    Internal use only.
    '''
    def __init__(self, prog, indent_increment=2, max_help_position=34, width=None):
        super().__init__(
            prog=prog,
            indent_increment=indent_increment,
            max_help_position=max_help_position,
            width=width)

    def _split_lines(self, text, width):
        # Allow explicit line breaks in help text using a ">>" sentinel.
        if text.startswith('>>'):
            return [line.strip() for line in text[2:].splitlines()]
        return super()._split_lines(text, width)

    def _format_action_invocation(self, action):
        # Strip hidden metavars so help alignment stays consistent.
        text = super()._format_action_invocation(action)
        return text.replace('\x08', '').rstrip()

    def _fill_text(self, text, width, indent):
        lines = functools.reduce(
            lambda x, y: x + y,
            map(
                lambda t: textwrap.wrap(
                    text=t.strip(),
                    width=width),
                text.split('\n\n')),
            [])

        for i in range(len(lines)):
            lines[i] = indent + lines[i]
            lines[i] = lines[i].split(' ')
            if lines[i][0] == '*' and lines[i][-1] == '+':
                lines[i] = '\r\n' + ' '.join(
                    lines[i][1:-1]) + '\n '
            elif lines[i][0] == '*':
                lines[i] = '\r\n' + ' '.join(
                    lines[i][1:])
            elif lines[i][-1] == '+':
                lines[i] = '\r' + ' '.join(
                    lines[i][:-1]) + '\n '
            else:
                lines[i] = '\r' + ' '.join(lines[i])
        return '\n'.join(lines)


def _print_header():
    """Print the CLI header banner once per process."""
    global header
    if not header:
        print()
        print(__doc__)
        print()
        header = True


def _print_footer():
    """Print the CLI footer banner once per process."""
    global footer
    if not footer:
        print(f'\nUTC {dt.datetime.utcnow()}\n')
        footer = True


def _parse_list_str(value):
    """Parse a comma-delimited string into a list of strings."""
    if value is None:
        return None
    if isinstance(value, list):
        return value
    if ',' in value:
        items = [v.strip() for v in value.split(',') if v.strip()]
        return items if items else None
    return value


def _parse_list_int(value):
    """Parse comma-delimited ints or a scalar int; pass through non-numeric strings."""
    if value is None:
        return None
    if isinstance(value, list):
        return value
    lowered = value.strip().lower()
    if lowered in ('none', 'null'):
        return None
    if ',' in value:
        items = [v.strip() for v in value.split(',') if v.strip()]
        if all(item.lstrip('-').isdigit() for item in items):
            return [int(item) for item in items]
        return value
    if value.lstrip('-').isdigit():
        return int(value)
    return value


def _dump_stats(stats, args):
    """Emit stats as JSON to stdout or to a file if requested."""
    if stats is None:
        return
    if args.stats_json:
        print(json.dumps(stats, indent=2, sort_keys=True, default=str))
    if args.stats_file:
        with open(args.stats_file, 'w', encoding='utf-8') as handle:
            json.dump(stats, handle, indent=2, sort_keys=True, default=str)


def _handle_result(result, args):
    """Extract stats from a module result and return a process exit code."""
    # Most module calls return (dataframe, stats); stats-only calls return a dict.
    stats = result[1] if isinstance(result, tuple) else result
    _dump_stats(stats, args)
    if isinstance(stats, dict) and 'status' in stats:
        return 0 if stats['status'] else 1
    return 0


def _add_common_options(parser, opt_group=None):
    """Register common CLI flags on the target parser/group."""
    target = opt_group if opt_group is not None else parser
    target.add_argument(
        '--stats-json',
        action='store_true',
        help='''>>[optional switch]
Print the stats dictionary as JSON to stdout.''')
    target.add_argument(
        '--stats-file',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Write the stats dictionary as JSON to this file path.''')
    target.add_argument(
        '--quiet',
        dest='verbose',
        action='store_false',
        help='''>>[optional switch]
Suppress verbose module output (sets verbose=False).''')
    parser.set_defaults(verbose=True)


def _add_manual(cmdpar):
    """Register the manual subcommand parser."""
    parser = cmdpar.add_parser(
        'manual',
        help='show module documentation',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    parser.add_argument(
        'topic',
        nargs='?',
        default=None,
        metavar='\b',
        help='''>>[optional string]
Module name for docs (e.g., barcode, primer). Use "topics" to list.''')
    return parser


def _print_manual(topic):
    """Print module or package documentation for manual command."""
    import oligopool as op

    def _strip_header(doc):
        if not doc:
            return doc
        lines = doc.strip().splitlines()
        if len(lines) >= 2 and lines[0].startswith('oligopool v') and lines[1].startswith('by '):
            lines = lines[2:]
            while lines and not lines[0].strip():
                lines = lines[1:]
        return '\n'.join(lines).strip()

    topics = {
        'background': background,
        'barcode': barcode,
        'primer': primer,
        'motif': motif,
        'spacer': spacer,
        'split': split,
        'pad': pad,
        'merge': merge,
        'revcomp': revcomp,
        'lenstat': lenstat,
        'final': final,
        'index': index,
        'pack': pack,
        'acount': acount,
        'xcount': xcount,
    }

    if topic is None:
        print(CLI_MANUAL.strip())
        return 0

    key = str(topic).strip().lower()
    if key in ('list', 'topics'):
        print('Available topics: {}'.format(', '.join(sorted(topics))))
        return 0
    if key in ('cli', 'manual'):
        print(CLI_MANUAL.strip())
        return 0
    if key in ('library', 'package'):
        doc = _strip_header(op.__doc__)
        if doc:
            print(doc)
            return 0
        print('No manual is available.')
        return 1

    target = topics.get(key)
    if target is None:
        print(f'No documentation available for "{topic}".')
        return 1

    doc = target.__doc__
    if doc:
        print(doc.strip())
        return 0

    print(f'No documentation available for "{topic}".')
    return 1


def _add_background(cmdpar):
    """Register the background subcommand parser."""
    parser = cmdpar.add_parser(
        'background',
        help='build background k-mer database',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to a CSV file with background sequences.
Required columns: ID, Sequence.''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length allowed between primers and background.
Range: 6 to 20.''')
    req.add_argument(
        '--output-directory',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output directory for the background k-mer database.
A ".oligopool.background" suffix is added if missing.''')
    _add_common_options(parser, opt)
    return parser


def _add_barcode(cmdpar):
    """Register the barcode subcommand parser."""
    parser = cmdpar.add_parser(
        'barcode',
        help='design constrained barcodes',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed oligo length for the design (>= 4).''')
    req.add_argument(
        '--barcode-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Length of the barcode to design (>= 4).''')
    req.add_argument(
        '--minimum-hamming-distance',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Minimum pairwise Hamming distance between barcodes (>= 1).''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length shared with existing oligos (>= 4).''')
    req.add_argument(
        '--barcode-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name to store barcodes (overwrites if present).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.barcode.csv" suffix is added if missing.''')
    opt.add_argument(
        '--barcode-type',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Barcode design mode: 0 = terminus optimized (fast), 1 = spectrum optimized (slow).
(default: 0)''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing left flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing right flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with ID and Exmotif columns.''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible barcode design.''')
    _add_common_options(parser, opt)
    return parser


def _add_primer(cmdpar):
    """Register the primer subcommand parser."""
    parser = cmdpar.add_parser(
        'primer',
        help='design constrained primers',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed oligo length for the design (>= 4).''')
    req.add_argument(
        '--primer-sequence-constraint',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
IUPAC degenerate sequence constraint for the primer.''')
    req.add_argument(
        '--primer-type',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Primer direction: 0 = forward, 1 = reverse.''')
    req.add_argument(
        '--minimum-melting-temperature',
        required=True,
        type=float,
        metavar='\b',
        help='''>>[required float]
Minimum primer melting temperature in degC (>= 25).''')
    req.add_argument(
        '--maximum-melting-temperature',
        required=True,
        type=float,
        metavar='\b',
        help='''>>[required float]
Maximum primer melting temperature in degC (<= 95).''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length shared with input oligos (>= 6).''')
    req.add_argument(
        '--primer-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name to store primers (overwrites if present).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.primer.csv" suffix is added if missing.''')
    opt.add_argument(
        '--paired-primer-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing the paired primer to match Tm against.''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing left flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing right flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with ID and Exmotif columns.''')
    opt.add_argument(
        '--background-directory',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Path to background k-mer database generated by background().''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible primer design.''')
    _add_common_options(parser, opt)
    return parser


def _add_motif(cmdpar):
    """Register the motif subcommand parser."""
    parser = cmdpar.add_parser(
        'motif',
        help='design or add motifs',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed oligo length for the design (>= 4).''')
    req.add_argument(
        '--motif-sequence-constraint',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
IUPAC degenerate sequence constraint or constant motif.''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length shared with input oligos (>= 4).''')
    req.add_argument(
        '--motif-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name to store motifs (overwrites if present).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.motif.csv" suffix is added if missing.''')
    opt.add_argument(
        '--motif-type',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Motif type: 0 = non-constant, 1 = constant.
(default: 0)''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing left flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing right flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with ID and Exmotif columns.''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible motif design.''')
    _add_common_options(parser, opt)
    return parser


def _add_spacer(cmdpar):
    """Register the spacer subcommand parser."""
    parser = cmdpar.add_parser(
        'spacer',
        help='design or insert spacers',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed oligo length for the design (>= 4).''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length shared with input oligos (>= 4).''')
    req.add_argument(
        '--spacer-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name to store spacers (overwrites if present).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.spacer.csv" suffix is added if missing.''')
    opt.add_argument(
        '--spacer-length',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Spacer length as an integer, a comma-separated list, or a CSV path
with ID and Length columns. If omitted, spacer length is auto-set
to reach the oligo length limit.''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing left flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing right flanking context for motif exclusion.
At least one of --left-context-column or --right-context-column is required.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with ID and Exmotif columns.''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible spacer design.''')
    _add_common_options(parser, opt)
    return parser


def _add_split(cmdpar):
    """Register the split subcommand parser."""
    parser = cmdpar.add_parser(
        'split',
        help='split oligos into fragments',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--split-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed length for split fragments (>= 4).''')
    req.add_argument(
        '--minimum-melting-temperature',
        required=True,
        type=float,
        metavar='\b',
        help='''>>[required float]
Minimum overlap melting temperature in degC (>= 4).''')
    req.add_argument(
        '--minimum-hamming-distance',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Minimum pairwise Hamming distance for overlaps (>= 1).''')
    req.add_argument(
        '--minimum-overlap-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Minimum overlap length in bp (>= 15).''')
    req.add_argument(
        '--maximum-overlap-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum overlap length in bp (>= 15).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.split.csv" suffix is added if missing.''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible split decisions.''')
    _add_common_options(parser, opt)
    return parser


def _add_pad(cmdpar):
    """Register the pad subcommand parser."""
    parser = cmdpar.add_parser(
        'pad',
        help='pad split oligos with primers',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed padded oligo length (>= 60).''')
    req.add_argument(
        '--split-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name containing split fragments to pad.''')
    req.add_argument(
        '--typeiis-system',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Type IIS enzyme name for pad excision (see pad() docs for options).''')
    req.add_argument(
        '--minimum-melting-temperature',
        required=True,
        type=float,
        metavar='\b',
        help='''>>[required float]
Minimum padding primer melting temperature in degC (>= 25).''')
    req.add_argument(
        '--maximum-melting-temperature',
        required=True,
        type=float,
        metavar='\b',
        help='''>>[required float]
Maximum padding primer melting temperature in degC (<= 95).''')
    req.add_argument(
        '--maximum-repeat-length',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum repeat length shared with input oligos (6 to 20).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.pad.csv" suffix is added if missing.''')
    opt.add_argument(
        '--random-seed',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Random seed for reproducible padding.''')
    _add_common_options(parser, opt)
    return parser


def _add_merge(cmdpar):
    """Register the merge subcommand parser."""
    parser = cmdpar.add_parser(
        'merge',
        help='merge elements into one column',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--merge-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name for merged DNA (overwrites if present).''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.merge.csv" suffix is added if missing.''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Left boundary column for merging; defaults to first column if omitted.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Right boundary column for merging; defaults to last column if omitted.''')
    _add_common_options(parser, opt)
    return parser


def _add_revcomp(cmdpar):
    """Register the revcomp subcommand parser."""
    parser = cmdpar.add_parser(
        'revcomp',
        help='reverse complement elements',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename. A ".oligopool.revcomp.csv" suffix is added if missing.''')
    opt.add_argument(
        '--left-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Left boundary column to reverse-complement; defaults to first column.''')
    opt.add_argument(
        '--right-context-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Right boundary column to reverse-complement; defaults to last column.''')
    _add_common_options(parser, opt)
    return parser


def _add_lenstat(cmdpar):
    """Register the lenstat subcommand parser."""
    parser = cmdpar.add_parser(
        'lenstat',
        help='compute length statistics',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--oligo-length-limit',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Maximum allowed oligo length for reporting (>= 4).''')
    _add_common_options(parser, opt)
    return parser


def _add_final(cmdpar):
    """Register the final subcommand parser."""
    parser = cmdpar.add_parser(
        'final',
        help='finalize library',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--input-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to input CSV with an ID column and DNA sequence columns.''')
    req.add_argument(
        '--output-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output CSV filename for the finalized library.
A ".oligopool.final.csv" suffix is added if missing.''')
    _add_common_options(parser, opt)
    return parser


def _add_index(cmdpar):
    """Register the index subcommand parser."""
    parser = cmdpar.add_parser(
        'index',
        help='index barcodes and associates',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--barcode-data',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to barcode CSV with an ID column and barcode column.''')
    req.add_argument(
        '--barcode-column',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Column name containing barcodes to index.''')
    req.add_argument(
        '--index-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output index filename (".oligopool.index" is appended if missing).''')
    opt.add_argument(
        '--barcode-prefix-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column with constant barcode prefix anchor.
At least one of --barcode-prefix-column or --barcode-suffix-column is required.''')
    opt.add_argument(
        '--barcode-suffix-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column with constant barcode suffix anchor.
At least one of --barcode-prefix-column or --barcode-suffix-column is required.''')
    opt.add_argument(
        '--barcode-prefix-gap',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Distance in bp between prefix anchor and barcode in reads.
(default: 0)''')
    opt.add_argument(
        '--barcode-suffix-gap',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Distance in bp between suffix anchor and barcode in reads.
(default: 0)''')
    opt.add_argument(
        '--associate-data',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Path to associate CSV with an ID column and associate column.''')
    opt.add_argument(
        '--associate-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column name containing associate sequences to index.''')
    opt.add_argument(
        '--associate-prefix-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column with constant associate prefix anchor.
At least one of --associate-prefix-column or --associate-suffix-column is required.''')
    opt.add_argument(
        '--associate-suffix-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column with constant associate suffix anchor.
At least one of --associate-prefix-column or --associate-suffix-column is required.''')
    _add_common_options(parser, opt)
    return parser


def _add_pack(cmdpar):
    """Register the pack subcommand parser."""
    parser = cmdpar.add_parser(
        'pack',
        help='pack fastq reads',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--r1-fastq-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Path to R1 FastQ file (gzipped ok).''')
    req.add_argument(
        '--r1-read-type',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
R1 orientation: 0 = forward, 1 = reverse.''')
    req.add_argument(
        '--pack-type',
        required=True,
        type=int,
        metavar='\b',
        help='''>>[required integer]
Pack storage mode: 0 = concatenated, 1 = merged.''')
    req.add_argument(
        '--pack-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output pack filename (".oligopool.pack" is appended if missing).''')
    opt.add_argument(
        '--minimum-r1-read-length',
        type=int,
        default=1,
        metavar='\b',
        help='''>>[optional integer]
Minimum R1 read length (default: 1).''')
    opt.add_argument(
        '--minimum-r1-read-quality',
        type=int,
        default=20,
        metavar='\b',
        help='''>>[optional integer]
Minimum average R1 read quality (default: 20).''')
    opt.add_argument(
        '--r2-fastq-file',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Path to R2 FastQ file (paired-end only).''')
    opt.add_argument(
        '--r2-read-type',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
R2 orientation: 0 = forward, 1 = reverse.''')
    opt.add_argument(
        '--minimum-r2-read-length',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Minimum R2 read length (default: None).''')
    opt.add_argument(
        '--minimum-r2-read-quality',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Minimum average R2 read quality (default: None).''')
    opt.add_argument(
        '--pack-size',
        type=float,
        default=3.0,
        metavar='\b',
        help='''>>[optional float]
Target million unique reads per pack (0.1 to 5.0).
(default: 3.0)''')
    opt.add_argument(
        '--core-count',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
CPU cores to use (0 = auto).''')
    opt.add_argument(
        '--memory-limit',
        type=float,
        default=0.0,
        metavar='\b',
        help='''>>[optional float]
GB of memory per core (0 = auto).''')
    _add_common_options(parser, opt)
    return parser


def _add_acount(cmdpar):
    """Register the acount subcommand parser."""
    parser = cmdpar.add_parser(
        'acount',
        help='association counting',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--index-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Index file path (".oligopool.index" is appended if missing).''')
    req.add_argument(
        '--pack-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Pack file path (".oligopool.pack" is appended if missing).''')
    req.add_argument(
        '--count-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output count matrix filename (".oligopool.acount.csv" is appended if missing).''')
    opt.add_argument(
        '--mapping-type',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Mapping mode: 0 = fast, 1 = sensitive.
(default: 0)''')
    opt.add_argument(
        '--barcode-errors',
        type=int,
        default=-1,
        metavar='\b',
        help='''>>[optional integer]
Max barcode errors (-1 = auto-infer).
(default: -1)''')
    opt.add_argument(
        '--associate-errors',
        type=int,
        default=-1,
        metavar='\b',
        help='''>>[optional integer]
Max associate errors (-1 = auto-infer).
(default: -1)''')
    opt.add_argument(
        '--core-count',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
CPU cores to use (0 = auto).''')
    opt.add_argument(
        '--memory-limit',
        type=float,
        default=0.0,
        metavar='\b',
        help='''>>[optional float]
GB of memory per core (0 = auto).''')
    _add_common_options(parser, opt)
    return parser


def _add_xcount(cmdpar):
    """Register the xcount subcommand parser."""
    parser = cmdpar.add_parser(
        'xcount',
        help='combinatorial counting',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--index-files',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Comma-separated list of index files for combinatorial counting.
Each path may omit the ".oligopool.index" suffix.''')
    req.add_argument(
        '--pack-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Pack file path (".oligopool.pack" is appended if missing).''')
    req.add_argument(
        '--count-file',
        required=True,
        type=str,
        metavar='\b',
        help='''>>[required string]
Output count matrix filename (".oligopool.xcount.csv" is appended if missing).''')
    opt.add_argument(
        '--mapping-type',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
Mapping mode: 0 = fast, 1 = sensitive.
(default: 0)''')
    opt.add_argument(
        '--barcode-errors',
        type=int,
        default=-1,
        metavar='\b',
        help='''>>[optional integer]
Max barcode errors (-1 = auto-infer).
(default: -1)''')
    opt.add_argument(
        '--associate-errors',
        type=int,
        default=-1,
        metavar='\b',
        help='''>>[optional integer]
Max associate errors (-1 = auto-infer).
(default: -1)''')
    opt.add_argument(
        '--core-count',
        type=int,
        default=0,
        metavar='\b',
        help='''>>[optional integer]
CPU cores to use (0 = auto).''')
    opt.add_argument(
        '--memory-limit',
        type=float,
        default=0.0,
        metavar='\b',
        help='''>>[optional float]
GB of memory per core (0 = auto).''')
    _add_common_options(parser, opt)
    return parser


def _get_parsers():
    """Create and return the top-level CLI parser."""
    mainpar = OligopoolParser(
        prog='oligopool',
        usage='oligopool COMMAND --argument=<value> ...',
        formatter_class=OligopoolFormatter,
        description=None,
        add_help=False,
        epilog='''Note: Run "oligopool COMMAND" to see command-specific options.''')
    mainpar._positionals.title = 'COMMANDS Available'
    mainpar._optionals.title = 'Optional Arguments'

    cmdpar = mainpar.add_subparsers(
        metavar=' ',
        dest='command',
        required=True)

    _add_manual(cmdpar)
    _add_background(cmdpar)
    _add_barcode(cmdpar)
    _add_primer(cmdpar)
    _add_motif(cmdpar)
    _add_spacer(cmdpar)
    _add_split(cmdpar)
    _add_pad(cmdpar)
    _add_merge(cmdpar)
    _add_revcomp(cmdpar)
    _add_lenstat(cmdpar)
    _add_final(cmdpar)
    _add_index(cmdpar)
    _add_pack(cmdpar)
    _add_acount(cmdpar)
    _add_xcount(cmdpar)

    return mainpar


def main(argv=None):
    """CLI entry point for Oligopool Calculator.

    Exit codes: 0 (success), 1 (runtime error or failed stats), 404 (argparse error).
    """
    arg_list = sys.argv[1:] if argv is None else list(argv)
    # Keep stdout clean for machine-readable JSON output.
    show_banner = '--stats-json' not in arg_list
    if show_banner:
        _print_header()
    parser = _get_parsers()
    args = parser.parse_args(arg_list)

    try:
        if args.command == 'manual':
            exit_code = _print_manual(args.topic)
            if show_banner and exit_code == 0:
                _print_footer()
            return exit_code
        if args.command == 'background':
            result = background(
                input_data=args.input_data,
                maximum_repeat_length=args.maximum_repeat_length,
                output_directory=args.output_directory,
                verbose=args.verbose)
        elif args.command == 'barcode':
            result = barcode(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                barcode_length=args.barcode_length,
                minimum_hamming_distance=args.minimum_hamming_distance,
                maximum_repeat_length=args.maximum_repeat_length,
                barcode_column=args.barcode_column,
                output_file=args.output_file,
                barcode_type=args.barcode_type,
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                excluded_motifs=_parse_list_str(args.excluded_motifs),
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'primer':
            result = primer(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                primer_sequence_constraint=args.primer_sequence_constraint,
                primer_type=args.primer_type,
                minimum_melting_temperature=args.minimum_melting_temperature,
                maximum_melting_temperature=args.maximum_melting_temperature,
                maximum_repeat_length=args.maximum_repeat_length,
                primer_column=args.primer_column,
                output_file=args.output_file,
                paired_primer_column=args.paired_primer_column,
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                excluded_motifs=_parse_list_str(args.excluded_motifs),
                background_directory=args.background_directory,
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'motif':
            result = motif(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                motif_sequence_constraint=args.motif_sequence_constraint,
                maximum_repeat_length=args.maximum_repeat_length,
                motif_column=args.motif_column,
                output_file=args.output_file,
                motif_type=args.motif_type,
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                excluded_motifs=_parse_list_str(args.excluded_motifs),
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'spacer':
            result = spacer(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                maximum_repeat_length=args.maximum_repeat_length,
                spacer_column=args.spacer_column,
                output_file=args.output_file,
                spacer_length=_parse_list_int(args.spacer_length),
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                excluded_motifs=_parse_list_str(args.excluded_motifs),
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'split':
            result = split(
                input_data=args.input_data,
                split_length_limit=args.split_length_limit,
                minimum_melting_temperature=args.minimum_melting_temperature,
                minimum_hamming_distance=args.minimum_hamming_distance,
                minimum_overlap_length=args.minimum_overlap_length,
                maximum_overlap_length=args.maximum_overlap_length,
                output_file=args.output_file,
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'pad':
            result = pad(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                split_column=args.split_column,
                typeIIS_system=args.typeiis_system,
                minimum_melting_temperature=args.minimum_melting_temperature,
                maximum_melting_temperature=args.maximum_melting_temperature,
                maximum_repeat_length=args.maximum_repeat_length,
                output_file=args.output_file,
                verbose=args.verbose,
                random_seed=args.random_seed)
        elif args.command == 'merge':
            result = merge(
                input_data=args.input_data,
                merge_column=args.merge_column,
                output_file=args.output_file,
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                verbose=args.verbose)
        elif args.command == 'revcomp':
            result = revcomp(
                input_data=args.input_data,
                output_file=args.output_file,
                left_context_column=args.left_context_column,
                right_context_column=args.right_context_column,
                verbose=args.verbose)
        elif args.command == 'lenstat':
            result = lenstat(
                input_data=args.input_data,
                oligo_length_limit=args.oligo_length_limit,
                verbose=args.verbose)
        elif args.command == 'final':
            result = final(
                input_data=args.input_data,
                output_file=args.output_file,
                verbose=args.verbose)
        elif args.command == 'index':
            result = index(
                barcode_data=args.barcode_data,
                barcode_column=args.barcode_column,
                index_file=args.index_file,
                barcode_prefix_column=args.barcode_prefix_column,
                barcode_suffix_column=args.barcode_suffix_column,
                barcode_prefix_gap=args.barcode_prefix_gap,
                barcode_suffix_gap=args.barcode_suffix_gap,
                associate_data=args.associate_data,
                associate_column=args.associate_column,
                associate_prefix_column=args.associate_prefix_column,
                associate_suffix_column=args.associate_suffix_column,
                verbose=args.verbose)
        elif args.command == 'pack':
            result = pack(
                r1_fastq_file=args.r1_fastq_file,
                r1_read_type=args.r1_read_type,
                pack_type=args.pack_type,
                pack_file=args.pack_file,
                minimum_r1_read_length=args.minimum_r1_read_length,
                minimum_r1_read_quality=args.minimum_r1_read_quality,
                r2_fastq_file=args.r2_fastq_file,
                r2_read_type=args.r2_read_type,
                minimum_r2_read_length=args.minimum_r2_read_length,
                minimum_r2_read_quality=args.minimum_r2_read_quality,
                pack_size=args.pack_size,
                core_count=args.core_count,
                memory_limit=args.memory_limit,
                verbose=args.verbose)
        elif args.command == 'acount':
            result = acount(
                index_file=args.index_file,
                pack_file=args.pack_file,
                count_file=args.count_file,
                mapping_type=args.mapping_type,
                barcode_errors=args.barcode_errors,
                associate_errors=args.associate_errors,
                callback=None,
                core_count=args.core_count,
                memory_limit=args.memory_limit,
                verbose=args.verbose)
        elif args.command == 'xcount':
            index_files = [v.strip() for v in args.index_files.split(',') if v.strip()]
            result = xcount(
                index_files=index_files,
                pack_file=args.pack_file,
                count_file=args.count_file,
                mapping_type=args.mapping_type,
                barcode_errors=args.barcode_errors,
                associate_errors=args.associate_errors,
                callback=None,
                core_count=args.core_count,
                memory_limit=args.memory_limit,
                verbose=args.verbose)
        else:
            parser.print_help()
            return 2
    except Exception as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 1

    exit_code = _handle_result(result, args)
    if show_banner and exit_code == 0:
        _print_footer()
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
