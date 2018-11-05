import logging
import pickle


def serialize(data, filename):
    with open(filename, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        # f.write(pickle.dumps(data))
        # f.write(jsonpickle.encode(data, keys=True))


def deserialize(filename):
    with open(filename, 'rb') as f:
        try:
            data = pickle.load(f)
            # data = jsonpickle.decode(f.read(), keys=True)
        except EOFError as e:
            logging.error("Deserialization error: couldn't unpicle file '{}'".format(filename))
            raise
    return data


def serialize_to_string(obj):
    import jsonpickle
    return jsonpickle.dumps(obj)


def deserialize_from_string(string):
    import jsonpickle
    return jsonpickle.loads(string)
