'''
Config file loader and validation for Oligopool CLI.

This module provides YAML config file support for:
1. Single-command execution from config files
2. Multi-step pipeline execution from a single config
3. CLI argument overrides for config values
4. Parallel/branching pipeline execution with DAG dependencies

Config precedence: CLI args > config file > defaults
'''

from collections import defaultdict

import yaml


def load_config(path):
    '''Load and parse a YAML config file.

    Parameters:
        path (str): Path to the YAML config file.

    Returns:
        dict: Parsed config dictionary.

    Raises:
        FileNotFoundError: If the config file doesn't exist.
        yaml.YAMLError: If the config file is invalid YAML.
    '''
    with open(path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    if config is None:
        return {}
    return config

def get_command_config(config, command):
    '''Extract the config section for a specific command.

    Parameters:
        config (dict): Full config dictionary.
        command (str): Command name (e.g., 'barcode', 'primer').

    Returns:
        dict: Config section for the command, or empty dict if not found.
    '''
    section = config.get(command, {})
    if section is None:
        return {}
    return dict(section)

def get_pipeline_steps(config):
    '''Extract ordered step list from pipeline config.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        list: List of step/command names in execution order.

    Example config:
        pipeline:
          name: "MPRA Design"
          steps:
            - primer
            - barcode
            - spacer
            - final
    '''
    pipeline = config.get('pipeline', {})
    if pipeline is None:
        return []
    steps = pipeline.get('steps', [])
    if steps is None:
        return []
    return list(steps)

def get_pipeline_name(config):
    '''Extract pipeline name from config.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        str: Pipeline name, or 'Unnamed Pipeline' if not specified.
    '''
    pipeline = config.get('pipeline', {})
    if pipeline is None:
        return 'Unnamed Pipeline'
    return pipeline.get('name', 'Unnamed Pipeline')

def validate_config_structure(config):
    '''Validate basic config structure and return warnings.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        list: List of warning messages (empty if valid).

    Notes:
        - This performs structural validation only.
        - Parameter validation is handled by the existing validation_parsing.py.
    '''
    warnings = []

    if not isinstance(config, dict):
        warnings.append('Config must be a YAML mapping/dictionary.')
        return warnings

    # Check pipeline section if present
    pipeline = config.get('pipeline')
    if pipeline is not None:
        if not isinstance(pipeline, dict):
            warnings.append("'pipeline' section must be a mapping/dictionary.")
        else:
            steps = pipeline.get('steps')
            if steps is not None and not isinstance(steps, list):
                warnings.append("'pipeline.steps' must be a list.")

    return warnings

def convert_config_keys_to_args(config_section):
    '''Convert config keys to CLI-compatible argument names.

    YAML config uses snake_case (e.g., 'barcode_length').
    CLI args use kebab-case (e.g., '--barcode-length').
    argparse converts kebab-case to snake_case attributes.

    This function ensures config keys match argparse attribute names.

    Parameters:
        config_section (dict): Config section for a command.

    Returns:
        dict: Config with keys converted to snake_case.
    '''
    converted = {}
    for key, value in config_section.items():
        # Convert any kebab-case keys to snake_case
        snake_key = key.replace('-', '_')
        converted[snake_key] = value
    return converted

def is_parallel_pipeline(config):
    '''Check if the pipeline uses parallel/DAG-style step definitions.

    Parallel pipelines use dict-style steps with name/command/after fields.
    Sequential pipelines use simple string lists.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        bool: True if pipeline uses parallel step definitions.

    Example parallel config:
        pipeline:
          steps:
            - name: fwd_primer
              command: primer
            - name: barcode
              command: barcode
              after: [fwd_primer]

    Example sequential config:
        pipeline:
          steps:
            - primer
            - barcode
    '''
    steps = get_pipeline_steps(config)
    if not steps:
        return False
    # If first step is a dict, it's parallel format
    return isinstance(steps[0], dict)

def parse_parallel_steps(config):
    '''Parse parallel pipeline steps into normalized step definitions.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        list: List of step dicts with keys: name, command, after, config_key

    Each returned step dict contains:
        - name: Step identifier (for dependencies)
        - command: The oligopool command to run
        - after: List of step names this step depends on
        - config_key: Key to look up step config (defaults to name)
    '''
    steps = get_pipeline_steps(config)
    parsed = []

    for step in steps:
        if isinstance(step, str):
            # Simple string format - treat as sequential
            parsed.append({
                'name': step,
                'command': step,
                'after': [],
                'config_key': step,
            })
        elif isinstance(step, dict):
            name = step.get('name')
            command = step.get('command', name)
            after = step.get('after', [])
            config_key = step.get('config', name)

            if name is None:
                raise ValueError(f"Pipeline step missing 'name': {step}")

            # Normalize 'after' to list
            if after is None:
                after = []
            elif isinstance(after, str):
                after = [after]

            parsed.append({
                'name': name,
                'command': command,
                'after': list(after),
                'config_key': config_key,
            })
        else:
            raise ValueError(f"Invalid pipeline step format: {step}")

    return parsed

def build_execution_levels(steps):
    '''Build execution levels from step dependencies (topological sort).

    Groups steps into levels where all steps in a level can run in parallel.
    Steps in level N+1 depend on steps in levels 0..N.

    Parameters:
        steps (list): Parsed step definitions from parse_parallel_steps().

    Returns:
        list: List of levels, where each level is a list of step dicts
              that can execute in parallel.

    Raises:
        ValueError: If there's a cycle in dependencies or unknown dependency.

    Example:
        Input steps:
          - {name: 'a', after: []}
          - {name: 'b', after: []}
          - {name: 'c', after: ['a', 'b']}
          - {name: 'd', after: ['c']}

        Output levels:
          [
            [{name: 'a', ...}, {name: 'b', ...}],  # Level 0: parallel
            [{name: 'c', ...}],                     # Level 1: after a,b
            [{name: 'd', ...}],                     # Level 2: after c
          ]
    '''
    if not steps:
        return []

    # Build lookup and validate dependencies
    step_by_name = {s['name']: s for s in steps}
    for step in steps:
        for dep in step['after']:
            if dep not in step_by_name:
                raise ValueError(
                    f"Step '{step['name']}' depends on unknown step '{dep}'"
                )

    # Kahn's algorithm for topological sort with level tracking
    in_degree = {s['name']: len(s['after']) for s in steps}
    dependents = defaultdict(list)
    for step in steps:
        for dep in step['after']:
            dependents[dep].append(step['name'])

    # Start with steps that have no dependencies
    levels = []
    current_level = [step_by_name[name] for name, deg in in_degree.items() if deg == 0]

    processed = set()
    while current_level:
        levels.append(current_level)
        next_level = []

        for step in current_level:
            processed.add(step['name'])
            for dependent_name in dependents[step['name']]:
                in_degree[dependent_name] -= 1
                if in_degree[dependent_name] == 0:
                    next_level.append(step_by_name[dependent_name])

        current_level = next_level

    # Check for cycles
    if len(processed) != len(steps):
        unprocessed = [s['name'] for s in steps if s['name'] not in processed]
        raise ValueError(
            f"Cycle detected in pipeline dependencies. "
            f"Steps involved: {', '.join(unprocessed)}"
        )

    return levels

def validate_parallel_pipeline(config):
    '''Validate a parallel pipeline config and return warnings/errors.

    Parameters:
        config (dict): Full config dictionary.

    Returns:
        list: List of warning/error messages (empty if valid).
    '''
    warnings = []

    try:
        steps = parse_parallel_steps(config)
    except ValueError as e:
        warnings.append(str(e))
        return warnings

    if not steps:
        warnings.append("No pipeline steps defined.")
        return warnings

    # Check for duplicate step names
    names = [s['name'] for s in steps]
    seen = set()
    for name in names:
        if name in seen:
            warnings.append(f"Duplicate step name: '{name}'")
        seen.add(name)

    # Check that config sections exist for each step
    for step in steps:
        config_key = step['config_key']
        if config_key not in config:
            warnings.append(
                f"Step '{step['name']}' references config section '{config_key}' "
                f"which is not defined."
            )

    # Try to build execution levels (checks for cycles)
    try:
        build_execution_levels(steps)
    except ValueError as e:
        warnings.append(str(e))

    return warnings
