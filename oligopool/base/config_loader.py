'''
YAML config parsing and validation helpers
for Oligopool CLI workflows.
Internal use only.
'''

import collections as cx

import yaml


# Parser and Setup Functions

def load_config(path):
    '''
    Load and parse YAML config file.
    Internal use only.

    :: path
       type - string
       desc - path to YAML config file
    '''
    with open(path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    if config is None:
        return {}
    return config

def get_command_config(config, command):
    '''
    Extract command config section.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
    :: command
       type - string
       desc - command name to lookup
    '''
    section = config.get(command, {})
    if section is None:
        return {}
    return dict(section)

def get_pipeline_steps(config):
    '''
    Extract pipeline step definitions.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
    '''
    pipeline = config.get('pipeline', {})
    if pipeline is None:
        return []
    steps = pipeline.get('steps', [])
    if steps is None:
        return []
    return list(steps)

def get_pipeline_name(config):
    '''
    Extract pipeline name string.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
    '''
    pipeline = config.get('pipeline', {})
    if pipeline is None:
        return 'Unnamed Pipeline'
    return pipeline.get('name', 'Unnamed Pipeline')

def validate_config_structure(config):
    '''
    Validate top-level config structure.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
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
    '''
    Normalize config keys for argparse.
    Internal use only.

    :: config_section
       type - dict
       desc - command config section
    '''
    converted = {}
    for key, value in config_section.items():
        # Convert any kebab-case keys to snake_case
        snake_key = key.replace('-', '_')
        converted[snake_key] = value
    return converted

def is_parallel_pipeline(config):
    '''
    Determine if pipeline is DAG-style.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
    '''
    steps = get_pipeline_steps(config)
    if not steps:
        return False
    # If first step is a dict, it's parallel format
    return isinstance(steps[0], dict)

def parse_parallel_steps(config):
    '''
    Parse and normalize pipeline steps.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
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
                raise ValueError(
                    "Pipeline step missing 'name': {}".format(step))

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
            raise ValueError(
                'Invalid pipeline step format: {}'.format(step))

    return parsed

def build_execution_levels(steps):
    '''
    Build execution levels from dependencies.
    Internal use only.

    :: steps
       type - list
       desc - normalized step definitions
    '''
    if not steps:
        return []

    # Build lookup and validate dependencies
    step_by_name = {s['name']: s for s in steps}
    for step in steps:
        for dep in step['after']:
            if dep not in step_by_name:
                raise ValueError(
                    "Step '{}' depends on unknown step '{}'".format(
                        step['name'], dep)
                )

    # Kahn's algorithm for topological sort with level tracking
    in_degree = {s['name']: len(s['after']) for s in steps}
    dependents = cx.defaultdict(list)
    for step in steps:
        for dep in step['after']:
            dependents[dep].append(step['name'])

    # Start with steps that have no dependencies
    levels = []
    current_level = [
        step_by_name[name] for name, deg in in_degree.items() \
            if deg == 0]

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
        unprocessed = [
            s['name'] for s in steps \
                if s['name'] not in processed]
        raise ValueError(
            'Cycle detected in pipeline dependencies. Steps involved: {}'.format(
                ', '.join(unprocessed))
        )

    return levels

def validate_parallel_pipeline(config):
    '''
    Validate DAG-style pipeline structure.
    Internal use only.

    :: config
       type - dict
       desc - full parsed config dictionary
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
            warnings.append("Duplicate step name: '{}'".format(name))
        seen.add(name)

    # Check that config sections exist for each step
    for step in steps:
        config_key = step['config_key']
        if config_key not in config:
            warnings.append(
                "Step '{}' references config section '{}' which is not defined.".format(
                    step['name'], config_key)
            )

    # Try to build execution levels (checks for cycles)
    try:
        build_execution_levels(steps)
    except ValueError as e:
        warnings.append(str(e))

    return warnings
