
import shutil


def print_header(text, indent=0):
    idt = "\t"*indent
    term_s = shutil.get_terminal_size((80, 20))
    n_free = max((term_s.columns - indent*8), 0)
    sep = "="*n_free
    print(idt + sep)
    print(idt + text)
    print(idt + sep)
