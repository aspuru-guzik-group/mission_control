INPUT_FILE_NAME = 'message.txt'
OUTPUT_FILE_NAME = 'echo.txt'

cfg_specs = {
    'EXAMPLE_ECHO_ECHO_CMD': {'default': 'echo'}
}
for key, spec in cfg_specs.items(): spec['output_key'] = key
