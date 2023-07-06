from scipy.fftpack import dct
import math
import numpy as np


def dct2(a):
    array = dct(dct(a.T, norm='forward').T, norm='forward')
    onedarray = []
    [onedarray.extend(x) for x in array[:]]
    return onedarray

def normL2(a):
    length = len(a)
    marray = a
    norm = 0
    for i in range(0, length):
        value = marray[i]
        norm += value * value

    norm = math.sqrt(norm)
    return norm


def normalizeNormL2(a):
    norm = normL2(a)
    marray = a
    if norm != 0:
        invnorm = 1 / norm
        length = len(a)

        for i in range(0, length):
            value = marray[i]
            marray[i] = value * invnorm

    return marray


def entropyShannonSubTriangle(array, width, xl, yl, xh, yh, pPerPixel):

    marray = array
    entropy = 0
    for y in range(yl, yh):

        yi = y * width

        xend = int(xh - y * xh / yh)
        for x in range(xl, xend):
            i = int(yi + x)
            value = marray[i]
            if value > 0:
                entropy += value * math.log(value)

            elif value < 0:
                entropy += -value * math.log(-value)

    entropy = -entropy

    if (pPerPixel):
        entropy = 2 * entropy / ((xh - xl) * (yh - yl))

    return entropy


def compute(pDoubleArrayImage, pPSFSupportDiameter):

    lWidth = np.shape(pDoubleArrayImage)[1]
    lHeight = np.shape(pDoubleArrayImage)[0]

    pDoubleArrayImage = dct2(pDoubleArrayImage)
    pDoubleArrayImage = normalizeNormL2(pDoubleArrayImage)

    lOTFSupportX = int(lWidth / pPSFSupportDiameter)
    lOTFSupportY = int(lHeight / pPSFSupportDiameter)
    lEntropy = entropyShannonSubTriangle(pDoubleArrayImage, lWidth, 0,
                                         0,
                                         lOTFSupportX,
                                         lOTFSupportY,
                                         True)

    return lEntropy
