

class WPM3(object):
    """ See http://web.udl.es/usuaris/q4374304/ """
    TAG = 'wpm3-co'

    def __init__(self, rundir):
        self.run_dir = rundir

    @staticmethod
    def command(input_filename):
        return ['WPM3-2015-co', input_filename]


class Maxino(object):
    TAG = 'maxino'

    def __init__(self, rundir):
        self.run_dir = rundir

    @staticmethod
    def command(input_filename):
        return ['maxino', input_filename]


class Openwbo(object):
    TAG = 'openwbo'

    def __init__(self, rundir):
        self.run_dir = rundir

    @staticmethod
    def command(input_filename):
        return ['open-wbo_static', input_filename]
