#simpleFSC, author John Heumann
import numpy as np
from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
import matplotlib.pyplot as plt

def fsc(vol1, vol2, fShells = False, edgeWidth = 2, winName=None, winParams=None, winCenter=None):
    if len(vol1.shape) != len(vol2.shape) or not np.all(vol1.shape == vol2.shape):
        raise ValueError("Both volumes must be the same size!")

    
    if  isinstance(fShells, bool):
        fShells = np.linspace(0, 0.5, int(np.floor(0.5 * np.max(np.array(vol1.shape)))))[1:]
    nX, nY, nZ = vol1.shape
    fX = ((np.arange(nX) - np.floor(nX / 2)) / nX)
    fY = ((np.arange(nY) - np.floor(nY / 2)) / nY)
    fZ = ((np.arange(nZ) - np.floor(nZ / 2)) / nZ)

    arrFX, arrFY, arrFZ = np.meshgrid(fY, fX, fZ)
    fMag = np.sqrt(arrFX**2 + arrFY**2 + arrFZ**2)

    nShells = len(fShells)
    fscc = np.zeros(nShells)
    nSamples = np.zeros(nShells)
    fLow = np.concatenate(([0], fShells[:-1]))

    if winName is not None:
        if winName.lower() == 'gaussian':
            win = gaussian_filter(np.ones(vol1.shape), winParams, mode='constant') #bad translation
            plt.imshow(win[20])
        elif winName.lower() == 'file':
            win = np.load(winParams)
        else:
            win_x = np.hanning(vol1.shape[0])
            win_y = np.hanning(vol1.shape[1])
            win_z = np.hanning(vol1.shape[2])
            win = win_x[:, None, None] * win_y[None, :, None] * win_z[None, None, :]
            
        edge_mask = np.ones_like(vol1)
        if edgeWidth > 0:
            edge_mask[edgeWidth:-edgeWidth+1, edgeWidth:-edgeWidth+1, edgeWidth:-edgeWidth+1] = 0
        edge_mean1 = vol1[edge_mask > 0.5].mean()
        edge_mean2 = vol2[edge_mask > 0.5].mean()
        vol1 = vol1 - edge_mean1
        vol2 = vol2 - edge_mean2
        vol1 = win * vol1
        vol2 = win * vol2
        vol1 += edge_mean1
        vol2 += edge_mean2

    VOL1 = fftshift(np.fft.fftn(vol1))
    VOL2 = fftshift(np.fft.fftn(vol2))

    for iShell in range(nShells):
        idxShell = (fMag >= fLow[iShell]) & (fMag < fShells[iShell])
        shell1 = VOL1[idxShell]
        shell2 = VOL2[idxShell]

        nrgShell1 = np.sum(np.abs(shell1)**2)
        nrgShell2 = np.sum(np.abs(shell2)**2)
        fscc[iShell] = np.real(np.dot(shell1, shell2.conj())) / np.sqrt(nrgShell1 * nrgShell2)
        nSamples[iShell] = np.sum(idxShell)

    return fscc, fShells
