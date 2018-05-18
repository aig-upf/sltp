
import shutil


def print_header(text, indent=0):
    term_s = shutil.get_terminal_size((80, 20))
    n_free = min(max((term_s.columns - indent*8), 0), 80)
    sep = "="*n_free
    print_lines([sep, text, sep], indent)


def print_lines(text, indent=0):
    text = [text] if isinstance(text, str) else text
    idt = "\t"*indent
    print(idt + ("\n"+idt).join(text))
