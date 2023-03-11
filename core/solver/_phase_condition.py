import numpy as np


def phase_condition(self):
    # parse and sort phase condition indices. range defined in json file is inclusive
    h_idx = []
    idx = self.prob.cont_params["continuation"]["phase_cond_index"]
    if idx and idx != "allvel":
        idx = idx.split(",")
        for i in range(len(idx)):
            if "-" in idx[i]:
                idxrange = idx[i].split("-")
                h_idx.extend(list(range(int(idxrange[0]), int(idxrange[1]) + 1)))
            else:
                h_idx.append(int(idx[i]))
        h_idx = sorted(set(h_idx))
    elif idx == "allvel":
        sizeX = len(self.X0)
        h_idx = list(range(sizeX//2, sizeX))    

    # create phase condition matrix h
    self.nphase = len(h_idx)
    self.h = np.zeros((self.nphase, len(self.X0)))
    self.h[list(range(self.nphase)), h_idx] = 1.0
