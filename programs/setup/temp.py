#! /usr/bin/env python3
import numpy as np

c = np.arange(24).reshape((4,6))
print(np.asarray(np.matrix(c)).flatten())
