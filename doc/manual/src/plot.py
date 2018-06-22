import numpy
import matplotlib

# matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": [
         r"\usepackage[utf8x]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{cmbright}",
         ],
    "font.family" : "serif",
    "font.size": 10,
})

import pylab

data = pylab.loadtxt("timings1", skiprows=1)
for i in range(2,5):
    data = pylab.minimum(data, pylab.loadtxt("timings"+str(i), skiprows=1))

pylab.figure(figsize=(5,2.8))
pylab.plot(data[:,0], data[:,1], 'o-k', label='operator()')
pylab.plot(data[:,0], data[:,2], 'o-b', label='evaluate()')
pylab.plot(data[:,0], data[:,3], 'o-r', label='std::function')
pylab.plot(data[:,0], data[:,4], 'o-g', label='virtual evaluate()')
pylab.ylim(ymax = 350, ymin = 0)
pylab.xlim(xmax = 16, xmin = 0)

pylab.legend()
# pylab.title("Title of Plot")
pylab.xlabel("range vector size $N$")
pylab.ylabel("ms")

pylab.tight_layout()

pylab.savefig("timings.pgf")
# pylab.show()
