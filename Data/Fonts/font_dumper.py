#!/usr/bin/env python
# Simple script to dump a TTF file into an array of unsigned chars that
# can be read into FreeType.

import argparse

parser = argparse.ArgumentParser(description='Dump TTF file to char array.')
parser.add_argument('--ttf-file', required=True, help='Name of TTF file.')
parser.add_argument('--output-file', required=True, help='Name of output file.')
parser.add_argument('--variable-name', default='raw_data', help='Name of variable for char array.'
                    '  Default=%(default)s.')
parser.add_argument('--string-name', default='ttf_font_data',
                    help='Name of variable for string array.'
                    '  Default=%(default)s.')
args = parser.parse_args()

with open(args.ttf_file, 'rb') as f:
  hexdata = f.read().hex()
num = 0

with open(args.output_file, 'w') as f:
  f.write('namespace {\n')
  f.write(f'const unsigned char {args.variable_name}[] = {{\n   ')
  for i in range(0, len(hexdata), 2):
    f.write(f' 0x{hexdata[i:i+2]},')
    num += 1
    if num == 12:
      f.write('\n   ')
      num = 0
  f.write('\n};\n')
  f.write('}  // namespace\n')
  f.write(f'const std::string {args.string_name}((const char *){args.variable_name},'
          f' (const char *){args.variable_name} + sizeof({args.variable_name}));')
