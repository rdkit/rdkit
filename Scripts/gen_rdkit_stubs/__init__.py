import sys
import os
import builtins
import importlib
import re
import pathlib
import subprocess
from multiprocessing.pool import ThreadPool
import tempfile
import shutil
import logging

logger = logging.getLogger(__name__)


WORKER_SCRIPT = "worker.py"
IMPORT_MODULES = re.compile(r"^\s*import\s+(.*)$")
FROM_IMPORT_MODULES = re.compile(r"^\s*from\s+(\S+)\s+import\s+(.*)$")
PYI_TO_PY = re.compile(r"\.pyi$")
RDKIT_MODULE_NAME = "rdkit"
RDBASE_MODULE_NAME = "rdBase"
INIT_PY = "__init__.py"

def rdkit_has_rdbase(rdkit_dir):
    """Returns True if rdkit_dir contains the rdBase Python module.

    Args:
        rdkit_dir (str): directory path

    Returns:
        bool: True if rdkit_dir contains the rdBase Python module
    """
    return any(f.startswith(RDBASE_MODULE_NAME) for f in os.listdir(rdkit_dir))

def purge_rdkit_source_dir_from_sys_path():
    """Remove the rdkit source dir from sys.path if present."""
    if os.path.isdir(RDKIT_MODULE_NAME) and not rdkit_has_rdbase(RDKIT_MODULE_NAME):
        abs_cwd = os.path.abspath(os.getcwd())
        indices_to_pop = sorted([i for i, p in enumerate(sys.path) if os.path.abspath(p) == abs_cwd], reverse=True)
        for i in indices_to_pop:
            sys.path.pop(i)

def find_rdkit_site_packages_path():
    site_packages_path = None
    for path_entry in sys.path:
        if not path_entry:
            continue
        rdkit_path = os.path.join(path_entry, RDKIT_MODULE_NAME)
        if os.path.isdir(rdkit_path) and rdkit_has_rdbase(rdkit_path):
            site_packages_path = path_entry
            break
    if site_packages_path is None:
        raise ValueError("Failed to find rdkit in PYTHONPATH")
    return site_packages_path

def find_rdkit_include_path():
    rdkit_include_path = None
    try:
        rdkit_include_path = find_rdkit_site_packages_path()
    except ValueError:
        pass
    if rdkit_include_path is not None:
        for i in range(3):
            rdkit_include_path = os.path.dirname(rdkit_include_path)
        rdkit_include_path = os.path.join(rdkit_include_path, "include")
        if not os.path.isdir(os.path.join(rdkit_include_path, RDKIT_MODULE_NAME)):
            rdkit_include_path = None
    return rdkit_include_path

def copy_stubs(src_entry, outer_dirs):
    """Copy src_entry to each directory in outer_dirs.
    If src_entry is a directory it will be recursively copied.
    The src_entry path is stripped off any leading directory
    including the "rdkit" directory.

    Args:
        src_entry (str): full path to a file or directory
        outer_dirs (list[str]): list of destination directory paths
    """
    base_src_entry = None
    for part in pathlib.Path(src_entry).parts:
        if base_src_entry is not None:
            base_src_entry.append(part)
        if part == RDKIT_MODULE_NAME:
            base_src_entry = []
    for outer_dir in outer_dirs:
        dst_entry = os.path.join(outer_dir, *base_src_entry)
        if os.path.isdir(src_entry):
            shutil.copytree(src_entry, dst_entry, dirs_exist_ok=True)
        elif os.path.isfile(src_entry):
            shutil.copyfile(src_entry, dst_entry)

def run_worker(module_name):
    """Worker function passed to multiprocessing.Pool.map.

    Args:
        module_name (str): name of the Python module we
        are working on

    Returns:
        tuple[str, str]: tuple with stdout and stderr.
        The stderr contents are only reported if return code is nonzero
    """
    out = ""
    err = ""
    args = [
        "--tempdir", run_worker.tempdir,
        "--module-name", module_name,
    ]
    if run_worker.keep_incorrect_staticmethods:
        args.append("--keep-incorrect-staticmethods")
    cmd = run_worker.cmd + args + ["--outer-dirs"] + run_worker.outer_dirs
    proc = subprocess.run(cmd, capture_output=True)
    if proc.returncode:
        msg = proc.stderr.decode("utf-8") or "(no error message)"
        cmd_as_str = " ".join(cmd)
        err = f"\"{cmd_as_str}\" failed with:\n{msg}"
    if proc.stdout:
        out = proc.stdout.decode("utf-8")
    return out, err

def parse_modules_to_set(modules):
    """Convert comma-separated list of Python module names to a set.

    Args:
        modules (str): comma-separated list of Python module names

    Returns:
        set: set of unique Python module names
    """
    return {m.strip() for m in modules.split(",")}

def concat_parent_child_module(parent_module, child_module):
    """Convert "from x import y" to "import x.y".
    If y is "*", "from x import *" becomes "import x".

    Args:
        parent_module (str): from module
        child_module (str): imported module

    Returns:
        str: dot-separated concatenated import
    """
    if child_module != "*":
        parent_module += "." + child_module
    return parent_module

def clear_stubs(outer_dir):
    """Remove all files and directories from "rdkit-stubs"
    except CMakeLists.txt.

    Args:
        outer_dir (_type_): _description_
    """
    for entry in os.listdir(outer_dir):
        if entry in (
            "CMakeLists.txt",
            "gen_rdkit_stubs.out",
            "gen_rdkit_stubs.err"
        ):
            continue
        entry = os.path.join(outer_dir, entry)
        if os.path.isdir(entry):
            shutil.rmtree(entry)
        else:
            os.remove(entry)

def generate_stubs_internal(modules, outer_dirs, args):
    concurrency = min(args.concurrency, len(modules))
    with tempfile.TemporaryDirectory() as tempdir:
        src_dir = os.path.join(tempdir, RDKIT_MODULE_NAME)
        run_worker.cmd = [sys.executable, os.path.join(os.path.dirname(__file__), WORKER_SCRIPT)]
        run_worker.tempdir = tempdir
        run_worker.keep_incorrect_staticmethods = args.keep_incorrect_staticmethods
        run_worker.outer_dirs = outer_dirs
        with ThreadPool(concurrency) as pool:
            res = pool.map(run_worker, modules)
        concat_out, concat_err = tuple(zip(*res))
        concat_out = "\n".join(out for out in concat_out if out)
        concat_err = "\n".join(err for err in concat_err if err)
        if concat_err:
            logger.critical(concat_err)
        if concat_out and args.verbose:
            logger.warning(concat_out)
        if os.path.isdir(src_dir):
            for f in os.listdir(src_dir):
                src_entry = os.path.join(src_dir, f)
                if os.path.exists(src_entry):
                    copy_stubs(src_entry, outer_dirs)

def generate_stubs(site_packages_path, args):
    """Generate RDKit stubs.

    Args:
        site_packages_path (str): full path from where RDKit modules are
        imported.
        output_dirs (list, optional): List of directories where a rdkit-stubs.
        directory is created. Defaults to [os.getcwd()].
        concurrency (int, optional): Number of CPUs used to generate stubs.
        Defaults to auto.
        verbose (bool, optional): Whether output should be verbose.
        Defaults to False.
    """
    output_dirs = args.output_dirs or [os.getcwd()]
    args.concurrency = args.concurrency or 1
    args.verbose = args.verbose or False
    args.keep_incorrect_staticmethods = args.keep_incorrect_staticmethods or False
    modules = {str(p.parent.relative_to(site_packages_path)).replace(os.sep, ".")
               for p in sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob(INIT_PY))}
    outer_dirs = []
    for output_dir in output_dirs:
        outer_dir = os.path.join(output_dir, f"{RDKIT_MODULE_NAME}-stubs")
        outer_dirs.append(outer_dir)
        if os.path.isdir(outer_dir):
            clear_stubs(outer_dir)
        elif os.path.isfile(outer_dir):
            os.remove(outer_dir)
        if not os.path.isdir(outer_dir):
            os.makedirs(outer_dir)
    generate_stubs_internal(modules, outer_dirs, args)
    pyi_path = pathlib.Path(outer_dirs[0])
    pyi_files = {PYI_TO_PY.sub(".py", str(p.relative_to(pyi_path))) for p in pyi_path.rglob("*.pyi")}
    modules = {str(p.relative_to(site_packages_path.joinpath(RDKIT_MODULE_NAME)))
               for p in (
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.pyd")) +
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.so")) +
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.py"))
    ) if p.name != INIT_PY}.difference(pyi_files)
    modules = {RDKIT_MODULE_NAME + "." + os.path.splitext(p.replace(os.sep, "."))[0] for p in modules}
    generate_stubs_internal(modules, outer_dirs, args)

class PythonParameter:
    """Class to store Python function signature parameters."""

    ARG1 = "arg1"
    ARG_DIGIT = re.compile(r"^arg\d+$")

    def __init__(self, arg_type, arg_name, arg_default=None):
        self.arg_type = arg_type
        self.arg_name = arg_name
        self.arg_default = arg_default

    def as_str(self):
        """Return this parameter as a .pyi signature string parameter.

        Returns:
            str: parameter: type as string, followed by optional
            default parameter
        """
        res = f"{self.arg_name}: {self.arg_type}"
        if self.arg_default:
            res += f" = {self.arg_default}"
        return res

    def is_arg1(self):
        """Return True if this parameter is an arg1 parameter.

        Returns:
            bool: True if this parameter is an arg1 parameter
        """
        return self.arg_name == self.ARG1

    def is_arg_digit(self):
        """Return True if this parameter is an arg# parameter.

        Returns:
            bool: True if this parameter is an arg# parameter
        """
        return self.ARG_DIGIT.match(self.arg_name) is not None

    def rename(self, new_name):
        """Rename this parameter to new_name.

        Args:
            new_name (str): new parameter name
        """
        self.arg_name = new_name


class ProcessDocLines:
    """A class to pre-process docstrings before feeding them to pybind11_stubgen."""

    PY_SIGNATURE_ARG_REGEX = re.compile(r"^\((\S+)\)([^\=]+)\=?(.*)?$")
    DEF_REGEX = re.compile(r"^([^(]+)(\s*\(\s*).*(\s*\)\s*->\s*)[^:]+:(\s*)$")
    PURGE_CPP_OBJECT_ANGLE_BRACKETS = re.compile(r"^(.*)<(\S*\.)?(\S+)\s*object\s*at\s*\S+\s*>(.*)$")
    PURGE_OPEN_SQUARE_BRACKET = re.compile(r"\[(?!\])")
    PURGE_CLOSE_SQUARE_BRACKET = re.compile(r"(?<!\[)\]")
    PROTECTIONS = {
        "[": "__OPEN_SQUARE_BRACKET_TAG__",
        "]": "__CLOSE_SQUARE_BRACKET_TAG__",
        "=": "__EQUALS_TAG__"
    }
    CPP_PYTHONIC_RETURN_TYPES = {
        "_listSt6vectorIiSaIiEE": "typing.Sequence[typing.Sequence[int]]",
        "_ROConformerSeq": "typing.Sequence[rdkit.Chem.Conformer]",
        "_ROQAtomSeq": "typing.Sequence[rdkit.Chem.QueryAtom]",
        "_vectd": "typing.Sequence[double]",
        "_vecti": "typing.Sequence[int]",
        "_vectj": "typing.Sequence[int]",
        "_vectN5RDKit13Abbreviations22AbbreviationDefinitionE": "rdkit.Chem.rdAbbreviations.AbbreviationDefinition]",
        "_vectN5RDKit9Chirality10StereoInfoE": "typing.Sequence[rdkit.Chem.StereoInfo]",
        "_vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE": "typing.Sequence[str]",
        "_vectSt6vectorIiSaIiEE": "typing.Sequence[typing.Sequence[int]]",
        "_vectSt6vectorIjSaIjEE": "typing.Sequence[typing.Sequence[int]]",
    }
    STD_MEMBER_FUNC_PARAM_NAMES = {
        "__add__": ("other",),
        "__and__": ("other",),
        "append": ("item",),
        "__call__": (),
        "__contains__": ("item",),
        "__copy__": (),
        "data": (),
        "__delitem__": ("item",),
        "__enter__": (),
        "__eq__": ("other",),
        "__exit__": ("exc_type", "exc_value", "traceback"),
        "extend": ("other",),
        "__getinitargs__": (),
        "__getitem__": ("item",),
        "__getstate__": (),
        "__iadd__": ("other",),
        "__iand__": ("other",),
        "__idiv__": ("other",),
        "__imul__": ("other",),
        "__init__": (),
        "__invert__": (),
        "__ior__": ("other",),
        "__isub__": ("other",),
        "__iter__": (),
        "key": (),
        "__len__": (),
        "__mul__": ("other",),
        "__ne__": ("other",),
        "__next__": (),
        "next": (),
        "__or__": ("other",),
        "__setitem__": ("item", "value"),
        "__setstate__": ("data",),
        "__str__": (),
        "__sub__": ("other",),
        "__truediv__": ("other",),
        "__xor__": ("other",),
    }
    OVERLOADED_FUNCTION_TAG = "Overloaded function."

    def __init__(self, module_name, keep_incorrect_staticmethods=False):
        self.module_name = module_name
        self.keep_incorrect_staticmethods = keep_incorrect_staticmethods
        self.num_overloads = 0
        self.overload_num = 0
        self.top_signature = None

    @classmethod
    def protect_quoted_square_brackets_and_equals(cls, arg):
        """Replaces [, ], = with string tags.

        Args:
            arg (str): raw single Python signature parameter

        Returns:
            str: single Python signature parameter with symbols
            replaced by string tags
        """
        open_quote = False
        protected_arg = ""
        for c in arg:
            if c == "'":
                open_quote ^= True
            elif open_quote:
                replacement = cls.PROTECTIONS.get(c, None)
                if replacement is not None:
                    c = replacement
            protected_arg += c
        return protected_arg

    @classmethod
    def deprotect_quoted_square_brackets_and_equals(cls, arg):
        """Restores [, ], = which were replaced by string tags.

        Args:
            arg (str): single Python signature parameter with symbols
            replaced by string tags

        Returns:
            str: single Python signature parameter with original symbols
        """
        for k, v in cls.PROTECTIONS.items():
            arg = arg.replace(v, k)
        return arg

    @classmethod
    def process_py_signature_arg(cls, arg):
        """Process single Python signature parameter.

        Args:
            arg (str): raw single Python signature parameter

        Returns:
            str: processed single Python signature parameter
        """
        arg = cls.PURGE_CPP_OBJECT_ANGLE_BRACKETS.sub(r"\1\3()\4", arg)
        arg = cls.protect_quoted_square_brackets_and_equals(arg)
        arg = cls.PURGE_OPEN_SQUARE_BRACKET.sub("", arg)
        arg = cls.PURGE_CLOSE_SQUARE_BRACKET.sub("", arg)
        arg = arg.replace(" ", "")
        py_signature_arg_match = cls.PY_SIGNATURE_ARG_REGEX.match(arg)
        assert py_signature_arg_match
        py_signature_arg_type = py_signature_arg_match.group(1)
        py_signature_arg_name = py_signature_arg_match.group(2)
        py_signature_arg_default = cls.deprotect_quoted_square_brackets_and_equals(py_signature_arg_match.group(3))
        return PythonParameter(py_signature_arg_type, py_signature_arg_name, py_signature_arg_default)

    @classmethod
    def find_def_match(cls, src_line):
        """Find Python function def in src_line.

        Args:
            src_line (str): single docstring line

        Returns:
            re.Match|None: the re.Match or None if no match
        """
        return cls.DEF_REGEX.match(src_line)

    def convert_to_valid_type(self, py_signature_ret):
        """Takes as input a return type which may be either a builtin
        or an RDKit class type. Builtins are returned unchanged, while
        RDKit class types are prefixed with the appropriate module such
        that an import can be added at the top of the .pyi file by
        pybind11_stubgen and the return type is correctly recognized
        by Pyright.

        Args:
            py_signature_ret (str): raw return type

        Returns:
            str: return type prefixed with the relevant Python module
            if needed
        """
        res = self.CPP_PYTHONIC_RETURN_TYPES.get(py_signature_ret, py_signature_ret)
        is_valid_type = hasattr(builtins, py_signature_ret)
        module_name_tmp = self.module_name
        while module_name_tmp and not is_valid_type:
            module_tmp = importlib.import_module(module_name_tmp)
            is_valid_type = hasattr(module_tmp, res)
            if is_valid_type:
                res = f"{module_name_tmp}.{py_signature_ret}"
            else:
                module_name_tmp = ".".join(module_name_tmp.split(".")[:-1])
        return res

    def correct_function_args(self, func_name, args):
        """Correct Python signatures.

        In spite of our efforts, boost::python will still generate
        incorrect docstrings for certain class methods, i.e., with
        "arg1" instead of the "self" parameter. This happens, among
        others, for those generated with boost::python::make_constructor.
        We cannot patch the C++ sources adding the "self" parameter
        as a boost::python::arg or the compiler will issue an error.
        Therefore, the only option is to post-process the generated
        .pyi files and replace "arg1" with "self". When we do so
        and the parameters which follow are "arg#", we also renumber
        them accordingly.

        Args:
            func_name (str): name of the function whose signature
            may needs to be corrected
            args (list[PythonParameter]): list of PythonParameter
            instances that may need to be corrected

        Returns:
            list[PythonParameter]: corrected list of PythonParameter
            instances
        """
        param_names = self.STD_MEMBER_FUNC_PARAM_NAMES.get(func_name, None)
        if param_names is not None and args and args[0].is_arg1():
            param_names_with_self = ["self", *param_names]
            param_idx = 0
            for i, arg in enumerate(args):
                param_name = None
                if i <= len(param_names):
                    param_name = param_names_with_self[i]
                elif arg.is_arg_digit():
                    param_idx += 1
                    param_name = f"arg{param_idx}"
                if param_name is not None:
                    arg.rename(param_name)
        return args

    def process_src_line(self, src_line):
        """Process single docstring line.
        Function def lines are processed, others are returned with
        no change.
        Currently, processing involves:
        * Processing function parameters to assign a type
          (extracted from the docstring itself)
        * Processing return type to prepend the Pythonmodule
          where the type is defined as needed
        * For overloaded functiona and methods, the first overload becomes
          the top-line signature, followed by a line with the
          OVERLOADED_FUNCTION_TAG. Then all the overloads follow as
          listed in the docstring prepended by a "#." tag where
          # is the 1-based overload number.
          This is the overload format that pybind11_stubgen expects
          and that will trigger the addition of the @typing.overload
          decorators as well as the per-overload docstring trimming

        Args:
            src_line (str): single docstring line

        Returns:
            str: processed docstring line
        """
        def_match = self.find_def_match(src_line)
        if def_match:
            func_name = def_match.group(1)
            overload_prefix = ""
            self.overload_num += 1
            if self.num_overloads > 1:
                overload_prefix = f"{self.overload_num}. "
            func_open_bracket = def_match.group(2)
            func_end_bracket_and_arrow = def_match.group(3)
            func_colon_to_end = def_match.group(4)
            py_signature_regex = re.compile(r"^\s*" + func_name + r"\s*\((.*)\)\s*->\s*(\S+)\s*:\s*$")
            py_signature_match = py_signature_regex.match(src_line)
            if py_signature_match:
                py_signature_args = py_signature_match.group(1)
                if py_signature_args:
                    py_signature_args = py_signature_args.split(", ")
                else:
                    py_signature_args = []
                py_signature_ret = py_signature_match.group(2)
                py_signature_ret = self.convert_to_valid_type(py_signature_ret)
                processed_args = [self.process_py_signature_arg(py_signature_arg) for py_signature_arg in py_signature_args]
                if not self.keep_incorrect_staticmethods:
                    processed_args = self.correct_function_args(func_name, processed_args)
                processed_args = ", ".join(arg.as_str() for arg in processed_args)
                src_line = f"{func_name}{func_open_bracket}{processed_args}{func_end_bracket_and_arrow}{py_signature_ret}{func_colon_to_end}"
                if self.top_signature is None:
                    self.top_signature = src_line
                if overload_prefix:
                    src_line = overload_prefix + src_line
        return src_line

    def process_doc_lines(self, doc_lines):
        """Process the raw docstring lines.
        * Count the number of overloads for the function described
          in the docstring
        * Trim any empty lines at the beginning of the docstring,
          as pybind11_stubgen does expects no empty lines
          at the top of the docstring or it will misbehave
        * If the function has >1 overload, prepend the first overload
          as top signature followed by a 2nd line bearing the
          OVERLOADED_FUNCTION_TAG string

        Args:
            doc_lines (list[str]): raw docstring lines

        Returns:
            list[str]: processed docstring lines
        """
        self.num_overloads = len(tuple(filter(None, [self.find_def_match(doc_line) for doc_line in doc_lines])))
        doc_lines = list(map(self.process_src_line, doc_lines))
        i = 0
        for i, doc_line in enumerate(doc_lines):
            if doc_line:
                break
        for _ in range(i):
            doc_lines.pop(0)
        if self.num_overloads > 1:
            doc_lines.insert(0, self.top_signature)
            doc_lines.insert(1, self.OVERLOADED_FUNCTION_TAG)
        return doc_lines

    @classmethod
    def process(cls, module_name, doc_lines, keep_incorrect_staticmethods=False):
        """Process the raw docstring lines.
        This is a convenience static function that creates an instance
        of ProcessDocLines and calls process_doc_lines() on it.

        Args:
            module_name (str): fully qualified Python module
            name the docstring belongs to
            doc_lines (list[str]): raw docstring lines
            keep_incorrect_staticmethods (bool): if true, incorrectly
            typed staticmethods are left unmodified

        Returns:
            list[str]: processed docstring lines
        """
        instance = cls(module_name, keep_incorrect_staticmethods)
        return instance.process_doc_lines(doc_lines)
