import logging
import os
import subprocess
import sys
import tempfile

BASEDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

tests = [
    ("gripper.py", "aaai_prob01_gc",
     ["total complexity: 8", "Numerical[free]", "Numerical[Not(Equal(at_g,at))]", "Numerical[Exists(carry,<universe>)]",
      "Numerical[Exists(at_g,at-robby)]"]),

    ("blocks.py", "aaai_clear_x_simple_hybrid",
     ["total complexity: 8",
      # "Boolean[holding]", (sometimes  Atom[handempty] is found instead)
      "Boolean[And(holding,Nominal(a))]",
      "Numerical[Exists(Star(on),Nominal(a))]"]),

    ("blocks.py", "aaai_clear_x_simple_hybrid_gc",
     ["total complexity: 8",
      # "Boolean[holding]", (sometimes  Atom[handempty] is found instead)
      "Boolean[And(clear_g,holding)]",
      "Numerical[Exists(Star(on),clear_g)]"]),

    ("blocks.py", "aaai_bw_on_x_y_completeness_opt",
    ["total complexity: 17",
     # "Boolean[holding]",
     "Boolean[And(holding,Nominal(a))]",
     "Numerical[Exists(Star(on),Nominal(b))]", "Numerical[Exists(Star(on),Nominal(a))]",
     "Boolean[And(Exists(on,Nominal(b)),Nominal(a))]"]),
]


def test(script, configuration, expected_output):
    cwd = os.path.join(BASEDIR, "experiments")
    with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as f:
        command = [os.path.join(cwd, script), configuration, "--all"]
        print('Calling "{}". Output redirected to "{}"'.format(' '.join(command), f.name))
        retcode = subprocess.call(command, stdout=f, stderr=f)
        if retcode:
            print("Experiment returned error code {}".format(retcode))
            sys.exit(-1)

        f.flush()
        f.seek(0)
        output = f.read()
        # print(output)
        for s in expected_output:
            if s not in output:
                print('Expected string "{}" not found in output. Check experiment output at "{}"'.format(s, f.name))
                sys.exit(-1)


def runtests():
    for script, configuration, expected_output in tests:
        test(script, configuration, expected_output)
    print("All tests OK")


if __name__ == "__main__":
    runtests()
