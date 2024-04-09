import numpy as np

def sliding_window(arr, window, step):
    out = []
    for x in range(window//2, len(arr)-window//2, step):
        if window%2 == 0:
            out.append(arr[x-window//2:x+window//2])
        else:
            out.append(arr[x-window//2:x+window//2+1])
    return np.array(out)
