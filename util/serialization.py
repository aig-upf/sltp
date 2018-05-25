
import jsonpickle
import pickle


def serialize(data, filename):
    with open(filename, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        # f.write(pickle.dumps(data))
        # f.write(jsonpickle.encode(data, keys=True))


def deserialize(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
        # data = jsonpickle.decode(f.read(), keys=True)
    return data
