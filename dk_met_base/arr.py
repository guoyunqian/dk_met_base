# _*_ coding: utf-8 _*_

"""
  Array manipulating functions.
"""

import numpy as np


def conform_dims(dims, r, ndim):
    """
    Expands an array or scalar so that it conforms to the shape of
      the given dimension sizes.

    :param dims: An array of dimension sizes of which r will be conformed to.
    :param r: An numpy array whose dimensions must be a subset of dims.
    :param ndim: An array of dimension indexes to indicate which dimension
                 sizes indicated by dims match the dimensions in r.
    :return: This function will create a new variable that has dimensions dims
             and the same type as r.
             The values of r will be copied to all of the other dimensions.

    :Example:
    >>> x = np.arange(12).reshape(3,4)
    >>> x = conform_dims([2,3,5,4,2],x,[1,3])
    >>> print(x.shape)
    (2, 3, 5, 4, 2)
    """

    # reshape r to conform the number of dimension
    sz = r.shape
    cdim = np.ones(len(dims), dtype=np.int)
    cdim[ndim] = sz
    rr = np.reshape(r, cdim)

    # repeat r to conform the dimension
    for i, item in enumerate(dims):
        if cdim[i] == 1:
            rr = np.repeat(rr, item, axis=i)

    # return
    return rr


def unshape(a):
    """
    Convert multiple dimension array to 2d array, keep 1st dim unchanged.

    :param a: array_like, > 2D.
    :return: 2d array, old shape

    >>> a = np.arange(40).reshape(2,4,5)
    >>> a.shape
     (2, 4, 5)
    >>> b, oldshape = tools.unshape(a)
    >>> b.shape
     (2, 20)
    >>> c = tools.deunshape(b, oldshape)
    >>> c.shape
     (2, 4, 5)

    """

    if np.ndim(a) < 2:
        raise ValueError("a must be at least 2 dimension")

    oldshape = a.shape
    array2d = a.reshape(oldshape[0], -1)
    return array2d, oldshape


def deunshape(a, oldshape):
    """
    restore a to old shape.

    :param a: array_like
    :param oldshape: return array shape
    :return: ndarray

    :example:
     >>> a = np.arange(40).reshape(2,4,5)
     >>> a.shape
     (2, 4, 5)
     >>> b, oldshape = unshape(a)
     >>> b.shape
     (2, 20)
     >>> c = deunshape(b, oldshape)
     >>> c.shape
     (2, 4, 5)
    """

    arraynd = a.reshape(oldshape)
    return arraynd


def expand(a, ndim, axis=0):
    """
    expand 1D array to ndim array.

    :param a: 1D array_like
    :param ndim: number of dimensions
    :param axis: position of 1D array
    :return: narray.

    :Example:
     >>> x = np.array([1, 2, 3])
     >>> y = expand(x, 3, axis=1)
     >>> y.shape
     (1, 3, 1)
     >>> y
     array([[[1],
             [2],
             [3]]])
    """

    if axis < 0:
        axis = ndim + axis
    res = np.asarray(a)
    if res.ndim != 1:
        raise ValueError("input array must be one dimensional array")
    idx = [x for x in range(ndim)]
    idx.remove(axis)
    for i in idx:
        res = np.expand_dims(res, axis=i)
    return res


def mrollaxis(a, axis, start=0):
    """
    numpy.rollaxis 's MaskedArray version.

    :param a: array_like
    :param axis: moved axis
    :param start: moved start position.
    :return: ndarray
    """

    if not hasattr(a, 'mask'):
        return np.rollaxis(a, axis, start=start)
    else:
        mask = np.ma.getmaskarray(a)
        data = np.ma.getdata(a)
        mask = np.rollaxis(mask, axis, start=start)
        data = np.rollaxis(data, axis, start=start)
        out = np.ma.asarray(data)
        out.mask = mask
        return out


def scale_vector(in_vector, min_range, max_range,
                 vector_min=None, vector_max=None):
    """
    This is a utility routine to scale the elements of
    a vector or an array into a given data range. nan values is not changed.

    :param in_vector: The input vector or array to be scaled.
    :param min_range: The minimum output value of the scaled vector.
    :param max_range: The maximum output value of the scaled vector.
    :param vector_min: Set this value to the minimum value of the vector,
                       before scaling (vector_min < vector).
                       The default value is Min(vector).
    :param vector_max: Set this value to the maximum value of the vector,
                       before scaling (vector_max < maxvalue).
                       The default value is Max(vector).
    :return: A vector or array of the same size as the input,
             scaled into the data range given by `min_range` and
             `max_range'. The input vector is confined to the data
             range set by `vector_min` and `vector_max` before
             scaling occurs.
    """

    # make sure numpy array
    vector = np.array(in_vector)

    # check keyword parameters
    if vector_min is None:
        vector_min = np.nanmin(vector)
    if vector_max is None:
        vector_max = np.nanmax(vector)

    # Calculate the scaling factors
    scale_factor = [(
        (min_range * vector_max) -
        (max_range * vector_min)) / (vector_max - vector_min),
        (max_range - min_range) / (vector_max - vector_min)]

    # return the scaled vector
    return vector * scale_factor[1] + scale_factor[0]


def matching(in_a, in_b, nan=True):
    """
    Keeping array a's values with b's sort.

    :param in_a: nd array.
    :param in_b: nd array.
    :param nan: do not involve nan values.
    :return: the same length a array.

    :Examples:
    >>> aa = np.array([3, 4, 2, 10, 7, 3, 6])
    >>> bb = np.array([ 5,  7,  3, 9, 6, np.nan, 11])
    >>> print(matching(aa, bb))
    """
    a = in_a.flatten()
    b = in_b.flatten()
    if nan:
        index = np.logical_and(np.isfinite(a), np.isfinite(b))
        a[index][np.argsort(b[index])] = np.sort(a[index])
    else:
        a[np.argsort(b)] = np.sort(a)
    a.shape = in_a.shape
    return a
