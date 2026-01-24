#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Oligopool Calculator command-line interface.

This module provides the `oligopool` and `op` entry points (see `pyproject.toml`).

Behavior notes:
  - `op` prints a compact command menu.
  - `op COMMAND` prints command-specific options (argparse errors show help).
  - `op complete` prints/installs argcomplete setup (banner-free; safe for `eval`/sourcing).
  - Banners are suppressed when `--stats-json` is used so stdout can be parsed as JSON.

Implementation notes:
  - `argcomplete.autocomplete(parser)` must run before any output is printed.
  - Exit codes: 0 (success), 1 (runtime error / failed stats), 404 (argparse error).
'''

import ast
import argparse
import datetime as dt
import difflib
import functools
import importlib
import json
import os
from pathlib import Path
import re
import sys
import textwrap

import argcomplete

from . import __version__


# Printed by `_print_header()`.
BANNER_TEXT = f'oligopool v{__version__}\nby ah'

MAIN_MENU_DESCRIPTION = (
    'Oligopool Calculator is a suite of algorithms for\n'
    'automated design and analysis of oligopool libraries.'
)

_DOC_SECTION_HEADERS = {
    'Required Parameters:',
    'Optional Parameters:',
    'Returns:',
    'Notes:',
}

COMMAND_TO_API = {
    'background': ('background', 'background'),
    'barcode': ('barcode', 'barcode'),
    'primer': ('primer', 'primer'),
    'motif': ('motif', 'motif'),
    'spacer': ('spacer', 'spacer'),
    'split': ('split', 'split'),
    'pad': ('pad', 'pad'),
    'merge': ('merge', 'merge'),
    'revcomp': ('revcomp', 'revcomp'),
    'verify': ('verify', 'verify'),
    'lenstat': ('lenstat', 'lenstat'),
    'final': ('final', 'final'),
    'index': ('index', 'index'),
    'pack': ('pack', 'pack'),
    'acount': ('acount', 'acount'),
    'xcount': ('xcount', 'xcount'),
}

_MODULE_AST_CACHE: dict[str, ast.AST] = {}
_FUNC_DOC_CACHE: dict[tuple[str, str], str] = {}

# In argcomplete mode, the CLI is executed on every tab press; avoid any work
# that is irrelevant to completion quality (e.g., docstring parsing for epilogs).
_ARGCOMPLETE_MODE = bool(
    os.environ.get('_ARGCOMPLETE')
    or ('COMP_LINE' in os.environ and 'COMP_POINT' in os.environ)
)


def _load_api_func(command: str):
    module_name, func_name = COMMAND_TO_API[command]
    module = importlib.import_module(f'.{module_name}', package=__package__)
    return getattr(module, func_name)


def _get_cached_module_ast(module_name: str) -> ast.AST:
    tree = _MODULE_AST_CACHE.get(module_name)
    if tree is not None:
        return tree
    module_path = Path(__file__).resolve().parent / f'{module_name}.py'
    source = module_path.read_text(encoding='utf-8')
    tree = ast.parse(source)
    _MODULE_AST_CACHE[module_name] = tree
    return tree


def _get_func_docstring(module_name: str, func_name: str) -> str:
    cached = _FUNC_DOC_CACHE.get((module_name, func_name))
    if cached is not None:
        return cached

    tree = _get_cached_module_ast(module_name)
    doc = ''
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            doc = ast.get_docstring(node, clean=True) or ''
            break
    _FUNC_DOC_CACHE[(module_name, func_name)] = doc
    return doc


def _extract_doc_notes(doc: str):
    if not doc:
        return []
    lines = doc.splitlines()

    start_idx = None
    for i, line in enumerate(lines):
        if line.strip() == 'Notes:':
            start_idx = i + 1
            break
    if start_idx is None:
        return []

    section_lines = []
    for line in lines[start_idx:]:
        if line.strip() in _DOC_SECTION_HEADERS:
            break
        section_lines.append(line.rstrip())

    bullets = []
    current = None
    for raw in section_lines:
        text = raw.strip()
        if not text:
            continue

        if text.startswith('- '):
            if current:
                bullets.append(current.strip())
            current = text
            continue

        if text.startswith('* '):
            sub = text[2:].strip()
            if current:
                current += f' {sub}'
            else:
                current = f'- {sub}'
            continue

        if current:
            current += f' {text}'
        else:
            current = text

    if current:
        bullets.append(current.strip())

    normalized = []
    for bullet in bullets:
        bullet = bullet.strip()
        # Some docstrings contain nested bullets like "- - foo"; strip all leading markers.
        bullet = re.sub(r'^(?:[*-]\s*)+', '', bullet).strip()
        bullet = re.sub(r'\s+', ' ', bullet).strip()
        bullet = re.sub(
            r'`([a-z][a-z0-9_]*)=([^`]+)`',
            lambda m: f"--{m.group(1).replace('_', '-')} {m.group(2).strip()}",
            bullet,
        )
        bullet = re.sub(
            r'`([a-z][a-z0-9_]*)`',
            lambda m: (
                f'--{m.group(1).replace("_", "-")}'
                if '_' in m.group(1)
                else m.group(0)
            ),
            bullet,
        )
        normalized.append(f'- {bullet}')
    return normalized


def _notes_epilog(command: str):
    if _ARGCOMPLETE_MODE:
        return None
    module_name, func_name = COMMAND_TO_API[command]
    notes = _extract_doc_notes(_get_func_docstring(module_name, func_name))
    if not notes:
        return None
    parts = ['Notes:']
    parts.extend(notes)
    return '\n\n'.join(parts)


CLI_MANUAL = '''
Oligopool CLI Manual

Usage:
  oligopool COMMAND --argument=<value> ...
  op COMMAND --argument=<value> ...
  oligopool manual [TOPIC]
  op manual [TOPIC]

Examples:
  op manual barcode
  op manual primer
  op manual topic
  op manual topics

Notes:
  - Topics include all CLI commands plus a few meta topics (run: op manual topics or op manual topic).
  - Design/transform commands require --output-file in CLI mode.
  - Run "oligopool COMMAND" to see command-specific options.
  - Use "oligopool manual topics" (or "oligopool manual topic") to list manual topics.
  - Use "oligopool cite" to print citation information.
  - Use "op complete --print-instructions" to enable tab-completion.
'''

CLI_COMPLETE_MANUAL = '''
Enable tab-completion (argcomplete)

One-time install (recommended):
  auto: op complete --install
  bash: op complete --install bash
   zsh: op complete --install zsh
  fish: op complete --install fish

Print a snippet instead:
  bash: eval "$(op complete --shell bash)"
   zsh: eval "$(op complete --shell zsh)"
  fish: op complete --shell fish | source
'''

CITATION_TEXT = '''
If you use Oligopool Calculator in your research, please cite:

Hossain A, Cetnar DP, LaFleur TL, McLellan JR, Salis HM.
Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements.
ACS Synth Biol. 2024;13(12):4218-4232. doi:10.1021/acssynbio.4c00661

Paper: https://pubs.acs.org/doi/10.1021/acssynbio.4c00661
'''

MANUAL_COMMAND_TOPICS = tuple(sorted(COMMAND_TO_API))

MANUAL_META_TOPICS = (
    'topics',
    'topic',
    'list',
    'cli',
    'manual',
    'library',
    'package',
    'complete',
    'completion',
)

MANUAL_TOPIC_CHOICES = tuple(sorted(set(MANUAL_COMMAND_TOPICS) | set(MANUAL_META_TOPICS)))

# Banner toggles to prevent duplicate prints within a single process.
_HEADER_PRINTED = False
_FOOTER_PRINTED = False


# === Argument parsing ===

class OligopoolParser(argparse.ArgumentParser):
    '''
    Custom Parser to show help messages on
    error. Internal use only.
    '''
    @staticmethod
    def _suggest(value, choices, limit=3):
        if not value or not choices:
            return None
        matches = difflib.get_close_matches(value, choices, n=limit, cutoff=0.6)
        return matches if matches else None

    def _get_subcommands(self):
        for action in self._actions:
            if isinstance(action, argparse._SubParsersAction):
                return sorted(action.choices.keys())
        return []

    def _get_option_strings(self):
        return sorted(self._option_string_actions.keys())

    def error(self, message):
        if message:
            message = message.strip()
            if not message.endswith(':'):
                print(f'error:\n{message}\n')

            if 'invalid choice:' in message:
                try:
                    bad = message.split('invalid choice:', 1)[1].strip()
                    bad = bad.split('(', 1)[0].strip().strip("'\"")
                except Exception:
                    bad = None
                suggestions = self._suggest(bad, self._get_subcommands())
                if suggestions:
                    print('did you mean: {}\n'.format(', '.join(suggestions)))

            if 'unrecognized arguments:' in message:
                try:
                    bad = message.split('unrecognized arguments:', 1)[1].strip()
                    bad = bad.split()[0].strip().strip("'\"")
                except Exception:
                    bad = None
                suggestions = self._suggest(bad, self._get_option_strings())
                if suggestions:
                    print('did you mean: {}\n'.format(', '.join(suggestions)))
        try:
            self.print_help()
            print()
        except Exception:
            pass
        sys.exit(404)

    def parse_known_args(self, args=None, namespace=None):
        '''Treat unknown args as a hard error for subparsers.'''
        args, argv = super().parse_known_args(args, namespace)
        if argv:
            argv_string = ' '.join(argv)
            if argv_string:
                self.error(f'unrecognized arguments: {argv_string}')
        return args, argv

    def format_help(self):
        '''Format help output with grouped command listing.'''
        text = super().format_help()
        lines = text.splitlines(True)
        try:
            header_idx = next(i for i, line in enumerate(lines) if line.startswith('COMMANDS Available:'))
        except StopIteration:
            return text

        first_cmd_idx = header_idx + 1
        while first_cmd_idx < len(lines) and lines[first_cmd_idx].strip() == '':
            first_cmd_idx += 1
        if first_cmd_idx >= len(lines) or not lines[first_cmd_idx].startswith('    '):
            return text

        cmd_idx = first_cmd_idx
        command_lines = []
        while cmd_idx < len(lines) and lines[cmd_idx].startswith('    '):
            if lines[cmd_idx].strip():
                command_lines.append(lines[cmd_idx])
            cmd_idx += 1
        while cmd_idx < len(lines) and lines[cmd_idx].strip() == '':
            cmd_idx += 1

        manual_line = next((line for line in command_lines if line.lstrip().startswith('manual')), None)
        cite_line = next((line for line in command_lines if line.lstrip().startswith('cite')), None)
        complete_line = next((line for line in command_lines if line.lstrip().startswith('complete')), None)
        if manual_line is None or cite_line is None or complete_line is None:
            return text

        middle_lines = [line for line in command_lines if line not in (manual_line, cite_line, complete_line)]
        gap_after = {'spacer', 'background', 'pad', 'revcomp', 'final'}
        middle_out = []
        for line in middle_lines:
            middle_out.append(line)
            try:
                cmd = line.strip().split()[0]
            except Exception:
                cmd = None
            if cmd in gap_after:
                middle_out.append('\n')

        out = []
        out.extend(lines[:header_idx + 1])
        out.append('\n')
        out.append(manual_line)
        out.append(cite_line)
        out.append('\n')
        out.extend(middle_out)
        out.append('\n')
        out.append(complete_line)
        out.append('\n')
        out.extend(lines[cmd_idx:])
        return ''.join(out)


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
        text = text.replace('\x08', '')
        # If a flag uses the hidden metavar sentinel (`\b`) with nargs, argparse may
        # render an extra " [ ...]" token and leave double-spaces after stripping.
        # Keep the invocation compact so the first help line stays on the same row.
        if getattr(action, 'metavar', None) == '\b':
            text = re.sub(r'\s*\[\s*\.\.\.\s*\]\s*$', '', text)
        text = re.sub(r'\s+', ' ', text).rstrip()
        return text

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


# === Banner helpers ===

def _print_header():
    '''Print the CLI header banner once per process.'''
    global _HEADER_PRINTED
    if not _HEADER_PRINTED:
        print()
        print(BANNER_TEXT)
        print()
        _HEADER_PRINTED = True


def _print_footer():
    '''Print the CLI footer banner once per process.'''
    global _FOOTER_PRINTED
    if not _FOOTER_PRINTED:
        # Use a fixed EST offset to keep the label stable (EST = UTC-05:00).
        now = dt.datetime.now(dt.timezone(dt.timedelta(hours=-5)))
        timestamp = now.strftime('%Y-%m-%d %I:%M:%S %p')
        print(f'\nEST {timestamp}\n')
        _FOOTER_PRINTED = True


# === Common parsing helpers ===

_SEQ_CONSTRAINT_MAX_LEN = 100000


def _eval_py_string_expr(expr, *, max_len=_SEQ_CONSTRAINT_MAX_LEN):
    '''Safely evaluate a small subset of Python string expressions.

    Allowed:
      - String literals: 'NNNN', "ATG"
      - Repetition: 'N'*20, 20*'N'
      - Concatenation: 'A'*10 + 'C'*10
      - Parentheses for grouping

    Disallowed: names, function calls, indexing, attributes, f-strings.
    '''
    try:
        tree = ast.parse(expr, mode='eval')
    except SyntaxError as exc:
        raise argparse.ArgumentTypeError(
            f'Invalid Python-style string expression: {expr!r}'
        ) from exc

    def _eval_node(node):
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (str, int)):
                return node.value
            raise argparse.ArgumentTypeError(
                f'Unsupported constant in expression: {expr!r}'
            )

        if isinstance(node, ast.BinOp):
            left = _eval_node(node.left)
            right = _eval_node(node.right)

            if isinstance(node.op, ast.Add):
                if not (isinstance(left, str) and isinstance(right, str)):
                    raise argparse.ArgumentTypeError(
                        f'Only string + string is supported: {expr!r}'
                    )
                combined_len = len(left) + len(right)
                if combined_len > max_len:
                    raise argparse.ArgumentTypeError(
                        f'Expanded sequence constraint is too long ({combined_len} > {max_len}).'
                    )
                return left + right

            if isinstance(node.op, ast.Mult):
                if isinstance(left, str) and isinstance(right, int):
                    seq, repeat = left, right
                elif isinstance(left, int) and isinstance(right, str):
                    seq, repeat = right, left
                else:
                    raise argparse.ArgumentTypeError(
                        f'Only string * int is supported: {expr!r}'
                    )
                if repeat < 0:
                    raise argparse.ArgumentTypeError(
                        f'Repetition count must be >= 0: {expr!r}'
                    )
                expanded_len = len(seq) * repeat
                if expanded_len > max_len:
                    raise argparse.ArgumentTypeError(
                        f'Expanded sequence constraint is too long ({expanded_len} > {max_len}).'
                    )
                return seq * repeat

            raise argparse.ArgumentTypeError(
                f'Unsupported operator in expression: {expr!r}'
            )

        raise argparse.ArgumentTypeError(
            f'Unsupported Python syntax in expression: {expr!r}'
        )

    out = _eval_node(tree.body)
    if not isinstance(out, str):
        raise argparse.ArgumentTypeError(
            f'Expression must evaluate to a string: {expr!r}'
        )
    if len(out) > max_len:
        raise argparse.ArgumentTypeError(
            f'Expanded sequence constraint is too long ({len(out)} > {max_len}).'
        )
    return out


def _parse_sequence_constraint(value):
    '''Parse an IUPAC constraint, expanding simple repeat/py-string forms.

    Examples (quote to avoid shell globbing):
      - N*20
      - GCC+N*20+CCG
      - 'N'*20
      - 'A'*10 + 'C'*10
    '''
    if value is None:
        return None
    value = value.strip()
    if not value:
        return value

    match = re.fullmatch(r'([A-Za-z]+)\s*\*\s*(\d+)', value)
    if match:
        seq, repeat_str = match.groups()
        repeat = int(repeat_str)
        if repeat < 0:
            raise argparse.ArgumentTypeError('Repetition count must be >= 0.')
        expanded_len = len(seq) * repeat
        if expanded_len > _SEQ_CONSTRAINT_MAX_LEN:
            raise argparse.ArgumentTypeError(
                f'Expanded sequence constraint is too long ({expanded_len} > {_SEQ_CONSTRAINT_MAX_LEN}).'
            )
        return seq * repeat

    # Shorthand concatenation form: GCC+N*20+CCG
    # (This is not Python syntax; it's a CLI convenience.)
    if ('+' in value) and ("'" not in value) and ('"' not in value):
        parts = [part.strip() for part in value.split('+')]
        if not parts or any(not part for part in parts):
            raise argparse.ArgumentTypeError(
                f'Invalid sequence constraint expression: {value!r}'
            )
        out = []
        total_len = 0
        for part in parts:
            part_match = re.fullmatch(r'([A-Za-z]+)\s*\*\s*(\d+)', part)
            if part_match:
                seq, repeat_str = part_match.groups()
                repeat = int(repeat_str)
                if repeat < 0:
                    raise argparse.ArgumentTypeError(
                        'Repetition count must be >= 0.'
                    )
                expanded = seq * repeat
            elif re.fullmatch(r'[A-Za-z]+', part):
                expanded = part
            else:
                raise argparse.ArgumentTypeError(
                    f'Invalid sequence constraint expression: {value!r}'
                )
            total_len += len(expanded)
            if total_len > _SEQ_CONSTRAINT_MAX_LEN:
                raise argparse.ArgumentTypeError(
                    f'Expanded sequence constraint is too long ({total_len} > {_SEQ_CONSTRAINT_MAX_LEN}).'
                )
            out.append(expanded)
        return ''.join(out)

    if ("'" in value) or ('"' in value):
        return _eval_py_string_expr(value)

    return value


def _parse_list_str(value):
    '''Parse a comma-delimited string into a list of strings.'''
    if value is None:
        return None
    if isinstance(value, list):
        out = []
        for item in value:
            if item is None:
                continue
            if not isinstance(item, str):
                out.append(item)
                continue
            if ',' in item:
                out.extend([v.strip() for v in item.split(',') if v.strip()])
            else:
                out.append(item.strip())
        out = [v for v in out if v]
        return out if out else None
    if ',' in value:
        items = [v.strip() for v in value.split(',') if v.strip()]
        return items if items else None
    return value.strip()


def _parse_list_int(value):
    '''Parse comma-delimited ints or a scalar int; pass through non-numeric strings.'''
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


def _parse_type_param(value):
    '''Parse type parameter accepting int or string.'''
    if value is None:
        return None
    value = str(value).strip()
    if value.lstrip('-').isdigit():
        return int(value)
    return value


def _looks_like_path(value):
    '''Return True if a CLI argument looks like a filesystem path.'''
    if value is None:
        return False
    value = str(value).strip()
    if not value:
        return False
    lowered = value.lower()
    if lowered.startswith(('./', '../', '~/', '~\\')):
        return True
    if any(sep in value for sep in ('/', '\\')):
        return True
    if lowered.endswith(('.csv', '.tsv')):
        return True
    return os.path.exists(value)


def _dump_stats(stats, args):
    '''Emit stats as JSON to stdout or to a file if requested.'''
    if stats is None:
        return
    if args.stats_json:
        print(json.dumps(stats, indent=2, sort_keys=True, default=str))
    if args.stats_file:
        with open(args.stats_file, 'w', encoding='utf-8') as handle:
            json.dump(stats, handle, indent=2, sort_keys=True, default=str)


def _handle_result(result, args):
    '''Extract stats from a module result and return a process exit code.'''
    # Most module calls return (dataframe, stats); stats-only calls return a dict.
    stats = result[1] if isinstance(result, tuple) else result
    _dump_stats(stats, args)
    if isinstance(stats, dict) and 'status' in stats:
        return 0 if stats['status'] else 1
    return 0


# === Subcommand parsing ===

def _add_common_options(parser, opt_group=None):
    '''Register common CLI flags on the target parser/group.'''
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
    '''Register the manual subcommand parser.'''
    parser = cmdpar.add_parser(
        'manual',
        help='show module documentation',
        description='Show module documentation from the Python API (like help(...)).',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    topic = parser.add_argument(
        'topic',
        nargs='?',
        default=None,
        metavar='\b',
        help='''>>[optional string]
Module name for docs (e.g., barcode, primer). Use "topics" to list.''')
    topic.completer = argcomplete.completers.ChoicesCompleter(MANUAL_TOPIC_CHOICES)
    return parser


def _add_cite(cmdpar):
    '''Register the cite subcommand parser.'''
    parser = cmdpar.add_parser(
        'cite',
        help='show citation information',
        description='Print citation information and the paper link.',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    return parser


def _add_complete(cmdpar):
    '''Register the complete subcommand parser.'''
    parser = cmdpar.add_parser(
        'complete',
        help='print or install shell completion',
        description='Print or install shell tab-completion setup (argcomplete).',
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    parser.add_argument(
        '--shell',
        choices=('bash', 'zsh', 'fish'),
        default=None,
        metavar='\b',
        help='''>>[optional string]
Target shell for completion output (bash, zsh, fish).''')
    parser.add_argument(
        '--install',
        nargs='?',
        const='auto',
        choices=('auto', 'bash', 'zsh', 'fish'),
        default=None,
        metavar='\b',
        help='''>>[optional string]
Install completion into your shell config (idempotent).
If omitted, auto-detect shell.''')
    parser.add_argument(
        '--print-instructions',
        action='store_true',
        help='''>>[optional switch]
Print short instructions instead of shell code.''')
    return parser


def _detect_shell():
    shell = os.environ.get('SHELL') or ''
    shell = os.path.basename(shell)
    if shell in ('bash', 'zsh', 'fish'):
        return shell
    return 'bash'


def _completion_snippet(shell, include_hint=False, lazy=False):
    hint_lines = []
    if include_hint:
        hint_lines = [
            '# Oligopool CLI tab-completion (argcomplete)',
            '# One-time install (recommended): op complete --install  (restart shell)',
            f'# Or: op complete --install {shell}',
            '# Show more help: op complete --print-instructions',
            '',
        ]
    if shell == 'bash':
        if lazy:
            return '\n'.join(hint_lines + [
                '__oligopool_argcomplete_init() {',
                '  [[ -n "${_OLIGOPOOL_ARGCOMPLETE_INSTALLED:-}" ]] && return 0',
                '  if command -v register-python-argcomplete >/dev/null 2>&1; then',
                '    eval "$(register-python-argcomplete op)"',
                '    eval "$(register-python-argcomplete oligopool)"',
                '    _OLIGOPOOL_ARGCOMPLETE_INSTALLED=1',
                '  fi',
                '}',
                '# Initialize now (if available), and retry on each prompt (e.g., after conda activate).',
                '__oligopool_argcomplete_init',
                'if [[ -n "${PROMPT_COMMAND:-}" ]]; then',
                '  PROMPT_COMMAND="__oligopool_argcomplete_init; ${PROMPT_COMMAND}"',
                'else',
                '  PROMPT_COMMAND="__oligopool_argcomplete_init"',
                'fi',
            ])
        return '\n'.join(hint_lines + [
            'if command -v register-python-argcomplete >/dev/null 2>&1; then',
            '  eval "$(register-python-argcomplete op)"',
            '  eval "$(register-python-argcomplete oligopool)"',
            'fi',
        ])
    if shell == 'zsh':
        if lazy:
            return '\n'.join(hint_lines + [
                'autoload -U bashcompinit && bashcompinit',
                '__oligopool_argcomplete_init() {',
                '  [[ -n "${_OLIGOPOOL_ARGCOMPLETE_INSTALLED:-}" ]] && return 0',
                '  if command -v register-python-argcomplete >/dev/null 2>&1; then',
                '    eval "$(register-python-argcomplete op)"',
                '    eval "$(register-python-argcomplete oligopool)"',
                '    _OLIGOPOOL_ARGCOMPLETE_INSTALLED=1',
                '  fi',
                '}',
                '# Initialize now (if available), and retry on each prompt (e.g., after conda activate).',
                '__oligopool_argcomplete_init',
                'precmd_functions+=(__oligopool_argcomplete_init)',
            ])
        return '\n'.join(hint_lines + [
            'autoload -U bashcompinit && bashcompinit',
            'if command -v register-python-argcomplete >/dev/null 2>&1; then',
            '  eval "$(register-python-argcomplete op)"',
            '  eval "$(register-python-argcomplete oligopool)"',
            'fi',
        ])
    if shell == 'fish':
        if lazy:
            return '\n'.join(hint_lines + [
                'function __oligopool_argcomplete_init --on-event fish_prompt',
                '  if set -q _OLIGOPOOL_ARGCOMPLETE_INSTALLED',
                '    return',
                '  end',
                '  if type -q register-python-argcomplete',
                '    register-python-argcomplete --shell fish op | source',
                '    register-python-argcomplete --shell fish oligopool | source',
                '    set -g _OLIGOPOOL_ARGCOMPLETE_INSTALLED 1',
                '    functions -e __oligopool_argcomplete_init',
                '  end',
                'end',
            ])
        return '\n'.join(hint_lines + [
            'register-python-argcomplete --shell fish op | source',
            'register-python-argcomplete --shell fish oligopool | source',
        ])
    raise ValueError(f'unsupported shell: {shell}')


def _upsert_rc_block(rc_path, shell):
    begin = '# >>> oligopool argcomplete >>>'
    end = '# <<< oligopool argcomplete <<<'
    snippet = _completion_snippet(shell, include_hint=False, lazy=True)
    block = f'{begin}\n{snippet}\n{end}\n'

    rc_path.parent.mkdir(parents=True, exist_ok=True)
    existing = rc_path.read_text(encoding='utf-8') if rc_path.exists() else ''
    if begin in existing and end in existing:
        pre = existing.split(begin, 1)[0]
        post = existing.split(end, 1)[1]
        updated = pre.rstrip('\n') + '\n\n' + block + '\n' + post.lstrip('\n')
    else:
        updated = existing.rstrip('\n') + '\n\n' + block
    rc_path.write_text(updated, encoding='utf-8')


def _install_completion(shell):
    home = Path.home()
    if shell == 'bash':
        _upsert_rc_block(home / '.bashrc', 'bash')
        return 0
    if shell == 'zsh':
        _upsert_rc_block(home / '.zshrc', 'zsh')
        return 0
    if shell == 'fish':
        comp_dir = home / '.config' / 'fish' / 'completions'
        comp_dir.mkdir(parents=True, exist_ok=True)
        snippet = _completion_snippet('fish', include_hint=False, lazy=True) + '\n'
        (comp_dir / 'op.fish').write_text(snippet, encoding='utf-8')
        (comp_dir / 'oligopool.fish').write_text(snippet, encoding='utf-8')
        return 0
    raise ValueError(f'unsupported shell: {shell}')


def _run_complete(shell, install, print_instructions):
    include_hint = shell is None
    shell = shell or _detect_shell()
    if print_instructions:
        print(CLI_COMPLETE_MANUAL.strip())
        return 0
    if install is not None:
        install_shell = shell if install == 'auto' else install
        _install_completion(install_shell)
        print(f'Installed completion for {install_shell}. Restart your shell.')
        return 0
    print(_completion_snippet(shell, include_hint=include_hint))
    return 0


def _print_manual(topic):
    '''Print module or package documentation for manual command.'''
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

    if topic is None:
        print(CLI_MANUAL.strip())
        return 0

    key = str(topic).strip().lower()
    if key in ('list', 'topic', 'topics'):
        command_topics = ', '.join(MANUAL_COMMAND_TOPICS)
        meta_topics = 'topic/topics/list, cli/manual, library/package, complete/completion'
        print('Available command topics:')
        print(textwrap.fill(command_topics, width=79, initial_indent='  ', subsequent_indent='  '))
        print('Available meta topics:')
        print(textwrap.fill(meta_topics, width=79, initial_indent='  ', subsequent_indent='  '))
        return 0
    if key in ('cli', 'manual'):
        print(CLI_MANUAL.strip())
        return 0
    if key in ('complete', 'completion'):
        print(CLI_COMPLETE_MANUAL.strip())
        return 0
    if key in ('library', 'package'):
        doc = _strip_header(op.__doc__)
        if doc:
            print(doc)
            return 0
        print('No manual is available.')
        return 1

    if key not in MANUAL_COMMAND_TOPICS:
        print(f'No documentation available for "{topic}".')
        return 1

    try:
        target = getattr(op, key)
    except Exception:
        print(f'No documentation available for "{topic}".')
        return 1

    doc = getattr(target, '__doc__', None)
    if doc:
        print(doc.strip())
        return 0

    print(f'No documentation available for "{topic}".')
    return 1


def _print_cite():
    '''Print citation information for the project.'''
    print(CITATION_TEXT.strip())
    return 0


# === Module subcommands ===

def _add_background(cmdpar):
    '''Register the background subcommand parser.'''
    parser = cmdpar.add_parser(
        'background',
        help='k-mer database for off-target screening',
        description='Build a background k-mer database for repeat exclusion.',
        epilog=_notes_epilog('background'),
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
    '''Register the barcode subcommand parser.'''
    parser = cmdpar.add_parser(
        'barcode',
        help='orthogonal barcodes with cross-set separation',
        description='Design constrained barcodes and write an output CSV.',
        epilog=_notes_epilog('barcode'),
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
        type=_parse_type_param,
        default=0,
        metavar='\b',
        help='''>>[optional int/string]
Barcode design mode: 0 or 'terminus' = terminus optimized (fast),
1 or 'spectrum' = spectrum optimized (slow). Aliases: 'term', 'fast', 'spec', 'slow'.
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
        '--patch-mode',
        dest='patch_mode',
        action='store_true',
        help='''>>[optional switch]
Fill missing values in an existing barcode column instead of creating a new column.''')
    opt.add_argument(
        '--incremental',
        dest='patch_mode',
        action='store_true',
        help=argparse.SUPPRESS)
    opt.add_argument(
        '--cross-barcode-columns',
        type=str,
        default=None,
        nargs='+',
        metavar='\b',
        help='''>>[optional string]
Comma-separated or space-separated column names whose barcodes must remain
distinct from the newly designed set. Must be used with --minimum-cross-distance.''')
    opt.add_argument(
        '--minimum-cross-distance',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Minimum pairwise Hamming distance between the new barcodes and the cross set.
Must be used with --cross-barcode-columns.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with an Exmotif column (ID column is allowed but not required).''')
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
    '''Register the primer subcommand parser.'''
    parser = cmdpar.add_parser(
        'primer',
        help='thermodynamic primers with Tm matching',
        description='Design constrained primers and write an output CSV.',
        epilog=_notes_epilog('primer'),
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
        type=_parse_sequence_constraint,
        metavar='\b',
        help='''>>[required string]
IUPAC degenerate sequence constraint for the primer.
Accepts Python-style expressions (quote as one argument), e.g. "'N'*20", or shorthand
concatenation like "GCC+N*20+CCG".''')
    req.add_argument(
        '--primer-type',
        required=True,
        type=_parse_type_param,
        metavar='\b',
        help='''>>[required int/string]
Primer direction: 0 or 'forward' for forward, 1 or 'reverse' for reverse.
Aliases: 'fwd', 'f', 'rev', 'r'.''')
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
        '--patch-mode',
        dest='patch_mode',
        action='store_true',
        help='''>>[optional switch]
Fill missing values in an existing primer column instead of creating a new column.''')
    opt.add_argument(
        '--incremental',
        dest='patch_mode',
        action='store_true',
        help=argparse.SUPPRESS)
    opt.add_argument(
        '--oligo-sets',
        type=str,
        default=None,
        nargs='+',
        metavar='\b',
        help='''>>[optional string]
Comma-separated or space-separated per-oligo set labels (one per input row), or a CSV path
with ID and OligoSet columns.
Primers are designed per set and checked for cross-set compatibility.''')
    opt.add_argument(
        '--paired-primer-column',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Column containing the paired primer to match Tm against.''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with an Exmotif column (ID column is allowed but not required).''')
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
    '''Register the motif subcommand parser.'''
    parser = cmdpar.add_parser(
        'motif',
        help='sequence motifs or constant anchors',
        description='Design or add motifs/anchors and write an output CSV.',
        epilog=_notes_epilog('motif'),
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
        type=_parse_sequence_constraint,
        metavar='\b',
        help='''>>[required string]
IUPAC degenerate sequence constraint or constant motif.
Accepts Python-style expressions (quote as one argument), e.g. "'N'*20", or shorthand
concatenation like "GCC+N*20+CCG".''')
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
        type=_parse_type_param,
        default=0,
        metavar='\b',
        help='''>>[optional int/string]
Motif type: 0 or 'per-variant' = non-constant, 1 or 'constant' = constant.
Aliases: 'var', 'non-constant', 'const', 'anchor'. (default: 0)''')
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
        '--patch-mode',
        dest='patch_mode',
        action='store_true',
        help='''>>[optional switch]
Fill missing values in an existing motif column instead of creating a new column.''')
    opt.add_argument(
        '--incremental',
        dest='patch_mode',
        action='store_true',
        help=argparse.SUPPRESS)
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with an Exmotif column (ID column is allowed but not required).''')
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
    '''Register the spacer subcommand parser.'''
    parser = cmdpar.add_parser(
        'spacer',
        help='neutral spacers to meet length targets',
        description='Design or insert spacers and write an output CSV.',
        epilog=_notes_epilog('spacer'),
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
        '--patch-mode',
        dest='patch_mode',
        action='store_true',
        help='''>>[optional switch]
Fill missing values in an existing spacer column instead of creating a new column.''')
    opt.add_argument(
        '--incremental',
        dest='patch_mode',
        action='store_true',
        help=argparse.SUPPRESS)
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with an Exmotif column (ID column is allowed but not required).''')
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
    '''Register the split subcommand parser.'''
    parser = cmdpar.add_parser(
        'split',
        help='break long oligos into overlapping fragments',
        description='Split long oligos into shorter ones and write an output CSV.',
        epilog=_notes_epilog('split'),
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
        '--separate-outputs',
        action='store_true',
        default=True,
        help='''>>[optional flag]
Write separate CSV files per split fragment (e.g., output.Split1.oligopool.split.csv).
Enabled by default in CLI mode. Use --no-separate-outputs to write a single combined file.''')
    opt.add_argument(
        '--no-separate-outputs',
        action='store_false',
        dest='separate_outputs',
        help='''>>[optional flag]
Write a single combined CSV with all Split columns instead of separate files.''')
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
    '''Register the pad subcommand parser.'''
    parser = cmdpar.add_parser(
        'pad',
        help='add excisable primer pads for scarless assembly',
        description='Add excisable primer pads to split fragments for scarless overlap-based assembly.',
        epilog=_notes_epilog('pad'),
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
    '''Register the merge subcommand parser.'''
    parser = cmdpar.add_parser(
        'merge',
        help='collapse contiguous columns',
        description='Merge oligo elements into one column and write an output CSV.',
        epilog=_notes_epilog('merge'),
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
    '''Register the revcomp subcommand parser.'''
    parser = cmdpar.add_parser(
        'revcomp',
        help='reverse complement a column range',
        description='Reverse complement spanning elements and write an output CSV.',
        epilog=_notes_epilog('revcomp'),
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


def _add_verify(cmdpar):
    '''Register the verify subcommand parser.'''
    parser = cmdpar.add_parser(
        'verify',
        help='QC constraints before synthesis',
        description='Verify a library DataFrame (length/motif/degeneracy checks; no output file).',
        epilog=_notes_epilog('verify'),
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
Path to input CSV with an ID column.''')
    opt.add_argument(
        '--oligo-length-limit',
        type=int,
        default=None,
        metavar='\b',
        help='''>>[optional integer]
Maximum allowed oligo length for reporting (>= 4).''')
    opt.add_argument(
        '--excluded-motifs',
        type=str,
        default=None,
        metavar='\b',
        help='''>>[optional string]
Comma-separated motifs or a CSV path with an Exmotif column (ID column is allowed but not required).''')
    _add_common_options(parser, opt)
    return parser


def _add_lenstat(cmdpar):
    '''Register the lenstat subcommand parser.'''
    parser = cmdpar.add_parser(
        'lenstat',
        help='length statistics and free-space check',
        description='Compute length statistics (prints results; no output file).',
        epilog=_notes_epilog('lenstat'),
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
    '''Register the final subcommand parser.'''
    parser = cmdpar.add_parser(
        'final',
        help='concatenate columns into synthesis-ready oligos',
        description='Finalize the library and write an output CSV.',
        epilog=_notes_epilog('final'),
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
    '''Register the index subcommand parser.'''
    parser = cmdpar.add_parser(
        'index',
        help='build barcode/associate index',
        description='Index barcodes (and optionally associates) into an index file.',
        epilog=_notes_epilog('index'),
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
    '''Register the pack subcommand parser.'''
    parser = cmdpar.add_parser(
        'pack',
        help='preprocess and deduplicate FastQ reads',
        description='Preprocess and pack FastQ reads for counting into a pack file.',
        epilog=_notes_epilog('pack'),
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
        type=_parse_type_param,
        metavar='\b',
        help='''>>[required int/string]
R1 orientation: 0 or 'forward' for forward, 1 or 'reverse' for reverse.
Aliases: 'fwd', 'f', 'rev', 'r'.''')
    req.add_argument(
        '--pack-type',
        required=True,
        type=_parse_type_param,
        metavar='\b',
        help='''>>[required int/string]
Pack storage mode: 0 or 'concatenated' for concatenated, 1 or 'merged' for merged.
Aliases: 'concatenate', 'concat', 'cat', 'joined', 'join', 'merge', 'assemble', 'assembled', 'asm'.''')
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
        type=_parse_type_param,
        default=None,
        metavar='\b',
        help='''>>[optional int/string]
R2 orientation: 0 or 'forward' for forward, 1 or 'reverse' for reverse.
Aliases: 'fwd', 'f', 'rev', 'r'.''')
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
    '''Register the acount subcommand parser.'''
    parser = cmdpar.add_parser(
        'acount',
        help='association counting (single index)',
        description='Execute association counting and write a count matrix CSV.',
        epilog=_notes_epilog('acount'),
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
        type=_parse_type_param,
        default=0,
        metavar='\b',
        help='''>>[optional int/string]
Mapping mode: 0 or 'fast' for fast, 1 or 'sensitive' for sensitive.
Aliases: 'quick', 'sens', 'accurate'. (default: 0)''')
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
    '''Register the xcount subcommand parser.'''
    parser = cmdpar.add_parser(
        'xcount',
        help='combinatorial counting (multiple indexes)',
        description='Execute combinatorial counting and write a count matrix CSV.',
        epilog=_notes_epilog('xcount'),
        usage=argparse.SUPPRESS,
        formatter_class=OligopoolFormatter,
        add_help=False)
    req = parser.add_argument_group('Required Arguments')
    opt = parser.add_argument_group('Optional Arguments')
    req.add_argument(
        '--index-files',
        required=True,
        type=str,
        nargs='+',
        metavar='\b',
        help='''>>[required string]
Comma-separated or space-separated list of index files for combinatorial counting.
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
        type=_parse_type_param,
        default=0,
        metavar='\b',
        help='''>>[optional int/string]
Mapping mode: 0 or 'fast' for fast, 1 or 'sensitive' for sensitive.
Aliases: 'quick', 'sens', 'accurate'. (default: 0)''')
    opt.add_argument(
        '--barcode-errors',
        type=int,
        default=-1,
        metavar='\b',
        help='''>>[optional integer]
Max barcode errors (-1 = auto-infer).
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
    '''Create and return the top-level CLI parser.'''
    mainpar = OligopoolParser(
        prog='oligopool',
        usage='oligopool COMMAND --argument=<value> ...',
        formatter_class=OligopoolFormatter,
        description=None,
        add_help=False,
        epilog='''Run "oligopool COMMAND" to see command-specific options.''')
    mainpar._positionals.title = 'COMMANDS Available'
    mainpar._optionals.title = 'Optional Arguments'

    cmdpar = mainpar.add_subparsers(
        metavar=' ',
        dest='command',
        required=True)

    _add_manual(cmdpar)
    _add_cite(cmdpar)
    _add_barcode(cmdpar)
    _add_primer(cmdpar)
    _add_motif(cmdpar)
    _add_spacer(cmdpar)
    _add_background(cmdpar)
    _add_split(cmdpar)
    _add_pad(cmdpar)
    _add_merge(cmdpar)
    _add_revcomp(cmdpar)
    _add_lenstat(cmdpar)
    _add_verify(cmdpar)
    _add_final(cmdpar)
    _add_index(cmdpar)
    _add_pack(cmdpar)
    _add_acount(cmdpar)
    _add_xcount(cmdpar)
    _add_complete(cmdpar)

    return mainpar


# === Entry point ===

def main(argv=None):
    '''CLI entry point for Oligopool Calculator.

    Exit codes: 0 (success), 1 (runtime error or failed stats), 404 (argparse error).
    '''
    arg_list = sys.argv[1:] if argv is None else list(argv)
    parser = _get_parsers()
    argcomplete.autocomplete(parser)
    # Keep stdout clean for machine-readable JSON output.
    command = arg_list[0] if arg_list else None
    show_banner = ('--stats-json' not in arg_list) and (command != 'complete')
    if show_banner:
        _print_header()
    if not arg_list:
        if show_banner:
            print(MAIN_MENU_DESCRIPTION)
            print()
        parser.print_help()
        if show_banner:
            _print_footer()
        return 0
    args = parser.parse_args(arg_list)

    try:
        match args.command:
            case 'manual':
                exit_code = _print_manual(args.topic)
                if show_banner:
                    _print_footer()
                return exit_code
            case 'cite':
                exit_code = _print_cite()
                if show_banner:
                    _print_footer()
                return exit_code
            case 'complete':
                exit_code = _run_complete(args.shell, args.install, args.print_instructions)
                if show_banner:
                    _print_footer()
                return exit_code
            case 'background':
                background = _load_api_func('background')
                result = background(
                    input_data=args.input_data,
                    maximum_repeat_length=args.maximum_repeat_length,
                    output_directory=args.output_directory,
                    verbose=args.verbose)
            case 'barcode':
                barcode = _load_api_func('barcode')
                cross_columns = _parse_list_str(args.cross_barcode_columns)
                if isinstance(cross_columns, str):
                    cross_columns = [cross_columns]
                if (cross_columns is None) ^ (args.minimum_cross_distance is None):
                    raise ValueError(
                        'cross_barcode_columns and minimum_cross_distance must be set together')
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
                    patch_mode=args.patch_mode,
                    cross_barcode_columns=cross_columns,
                    minimum_cross_distance=args.minimum_cross_distance,
                    excluded_motifs=_parse_list_str(args.excluded_motifs),
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'primer':
                primer = _load_api_func('primer')
                oligo_sets = args.oligo_sets
                if isinstance(oligo_sets, list):
                    oligo_sets = _parse_list_str(oligo_sets)
                    if isinstance(oligo_sets, list) and len(oligo_sets) == 1 and _looks_like_path(oligo_sets[0]):
                        oligo_sets = oligo_sets[0]
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
                    left_context_column=args.left_context_column,
                    right_context_column=args.right_context_column,
                    patch_mode=args.patch_mode,
                    oligo_sets=oligo_sets,
                    paired_primer_column=args.paired_primer_column,
                    excluded_motifs=_parse_list_str(args.excluded_motifs),
                    background_directory=args.background_directory,
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'motif':
                motif = _load_api_func('motif')
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
                    patch_mode=args.patch_mode,
                    excluded_motifs=_parse_list_str(args.excluded_motifs),
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'spacer':
                spacer = _load_api_func('spacer')
                result = spacer(
                    input_data=args.input_data,
                    oligo_length_limit=args.oligo_length_limit,
                    maximum_repeat_length=args.maximum_repeat_length,
                    spacer_column=args.spacer_column,
                    output_file=args.output_file,
                    spacer_length=_parse_list_int(args.spacer_length),
                    left_context_column=args.left_context_column,
                    right_context_column=args.right_context_column,
                    patch_mode=args.patch_mode,
                    excluded_motifs=_parse_list_str(args.excluded_motifs),
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'split':
                split = _load_api_func('split')
                result = split(
                    input_data=args.input_data,
                    split_length_limit=args.split_length_limit,
                    minimum_melting_temperature=args.minimum_melting_temperature,
                    minimum_hamming_distance=args.minimum_hamming_distance,
                    minimum_overlap_length=args.minimum_overlap_length,
                    maximum_overlap_length=args.maximum_overlap_length,
                    output_file=args.output_file,
                    separate_outputs=args.separate_outputs,
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'pad':
                pad = _load_api_func('pad')
                result = pad(
                    input_data=args.input_data,
                    oligo_length_limit=args.oligo_length_limit,
                    split_column=args.split_column,
                    typeIIS_system=args.typeiis_system,
                    minimum_melting_temperature=args.minimum_melting_temperature,
                    maximum_melting_temperature=args.maximum_melting_temperature,
                    maximum_repeat_length=args.maximum_repeat_length,
                    output_file=args.output_file,
                    random_seed=args.random_seed,
                    verbose=args.verbose)
            case 'merge':
                merge = _load_api_func('merge')
                result = merge(
                    input_data=args.input_data,
                    merge_column=args.merge_column,
                    output_file=args.output_file,
                    left_context_column=args.left_context_column,
                    right_context_column=args.right_context_column,
                    verbose=args.verbose)
            case 'revcomp':
                revcomp = _load_api_func('revcomp')
                result = revcomp(
                    input_data=args.input_data,
                    output_file=args.output_file,
                    left_context_column=args.left_context_column,
                    right_context_column=args.right_context_column,
                    verbose=args.verbose)
            case 'verify':
                verify = _load_api_func('verify')
                result = verify(
                    input_data=args.input_data,
                    oligo_length_limit=args.oligo_length_limit,
                    excluded_motifs=_parse_list_str(args.excluded_motifs),
                    verbose=args.verbose)
            case 'lenstat':
                lenstat = _load_api_func('lenstat')
                result = lenstat(
                    input_data=args.input_data,
                    oligo_length_limit=args.oligo_length_limit,
                    verbose=args.verbose)
            case 'final':
                final = _load_api_func('final')
                result = final(
                    input_data=args.input_data,
                    output_file=args.output_file,
                    verbose=args.verbose)
            case 'index':
                index = _load_api_func('index')
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
            case 'pack':
                pack = _load_api_func('pack')
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
            case 'acount':
                acount = _load_api_func('acount')
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
            case 'xcount':
                xcount = _load_api_func('xcount')
                result = xcount(
                    index_files=_parse_list_str(args.index_files),
                    pack_file=args.pack_file,
                    count_file=args.count_file,
                    mapping_type=args.mapping_type,
                    barcode_errors=args.barcode_errors,
                    callback=None,
                    core_count=args.core_count,
                    memory_limit=args.memory_limit,
                    verbose=args.verbose)
            case _:
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
