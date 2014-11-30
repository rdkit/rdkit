#! make_templates.py
"""Parse a C++ input and try to generate SWIG %template declarations so that this
    template <typename T1, typename T2>
    DoubleVect
    OnBitProjSimilarity(const T1& bv1,const T2& bv2);
becomes this
    %template(OnBitProjSimilarityEBV) OnBitProjSimilarity<ExplicitBitVect,ExplicitBitVect>;
"""
from __future__ import print_function

import re
template_match = re.compile(r"""template\s*\<(.+)\>\s*.*\s+(\w+)\s*\(.*""")

def make_entry(name, suffix,  template_types):
    template_part = '<' + ','.join(template_types) + '>'
    return '%template(' + name + suffix + ') ' + name + template_part + ';'
    
def handle_template_line(template_line,  suffix,  template_type):
    # print template_line
    match = template_match.search(template_line)
    types = match.group(1)
    name = match.group(2)
    template_types = len(types.split(',')) * [template_type]
    print(make_entry(name,  suffix,  template_types))
    
def find_declarations(file_name,  suffix,  template_type):
    f = open(file_name)
    in_template = False
    template_line = ''
    for line in f:
        if in_template:
            template_line += ' ' + line.strip()
            if line.find(')') >= 0:
                handle_template_line(template_line,  suffix,  template_type)
                in_template = False
        elif line.strip().startswith('template'):
            in_template = True
            template_line = line.strip()
    f.close()
    
