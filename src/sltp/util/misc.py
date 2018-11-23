
def update_dict(d, **kwargs):
    res = d.copy()
    res.update(kwargs)
    return res


def project_dict(obj, idxs):
    """ Project a given dictionary/set/list to only a subset of its elements, given by idxs """
    if isinstance(obj, dict):
        return {k: obj[k] for k in idxs}
    elif isinstance(obj, set):
        return {k for k in idxs if k in obj}
    elif isinstance(obj, set):
        return [o for o in obj if o in idxs]
    raise RuntimeError("Cannot project unexpected data structure '{}'".format(type(obj)))
