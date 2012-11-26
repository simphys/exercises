    xfolded = x.copy()
    xfolded -= np.floor(x/L)*L  #BECAUSE OF PBC
