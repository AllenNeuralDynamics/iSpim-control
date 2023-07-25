from scipy.fftpack import dct
from math import sqrt
import numpy as np


def dct2(a):
    array = dct(dct(a.T, norm='forward').T, norm='forward')
    onedarray = array.flatten()
    return onedarray

def normL2(a):
    a = np.sum(np.square(a))
    norm = sqrt(a)
    return norm


def normalizeNormL2(a):
    norm = normL2(a)
    if norm != 0:
        a = a * (1 / norm)
    return a


def entropyShannonSubTriangle(a, width, xl, yl, xh, yh, pPerPixel):
    a = abs(a)
    entropy = 0
    for y in range(yl, yh):
        yi = y * width
        xend = int(xh - y * xh / yh)
        entropy -= np.sum(a[xl + yi:xend + yi] * np.log(a[xl + yi:xend + yi]))

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
