Code for Synthesizing Energy-Bounded Cages
==========================================

Energy-bounded cages (EBCs) are a robust approach to robot grasping and pushing in the presence of forces such as gravity and friction. This repository contains the research code to analyze and synthesize EBCs.

Publications
-------------
This code accompanies the following publications.

[1] Synthesis of Energy-Bounded Planar Caging Grasps using Persistent Homology Jeffrey Mahler, Florian Pokorny, Sherdil Niyaz, Ken Goldberg. Workshop on the Algorithmic Foundations of Robotics, 2016.

[2] Energy-Bounded Caging: Formal Definition and 2D Energy Lower Bound Algorithm Based on Weighted Alpha Shapes Jeffrey Mahler, Florian Pokorny, Zoe McCarthy, A Frank van der Stappen, Ken Goldberg. IEEE Robotics and Automation Letters, 2016.  

Please cite [1] if you use this code in a research publication.

Installation
------------
The code is known to run on Ubuntu 14.04 LTS.

The following libraries are required for installation and usage:
* CGAL
* SOLID collision detection
* QHULL
* CGAL's geomview stream (will only work in linux)

Install the dependencies above as well as 3rdparty/libccd, then build using the standard cmake workflow:

```mkdir build && cd build
cmake ..
make -j4
```

If successful, the binary will be at build/src/alpha.

Data
----
Example data can be downloaded [here](https://github.com/BerkeleyAutomation/caging/raw/gh-pages/data.zip). To request additional data, email jmahler@berkeley.edu with the subject line "Request for Energy-Bounded Caging Data."

Authors
-------
Jeffrey Mahler: jmahler@berkeley.edu

Sherdil Niyaz: sniyaz@berkeley.edu

Zoe McCarthy: zoemccarthy12@gmail.com

AUTOLAB, UC Berkeley


Note
----
This project is ongoing research on using alpha shapes and persistent homology for proving
path existence / nonexistence for the anaysis and synthetis of robot caging configurations.
It was adapted from Zoe McCarthy's work at UIUC on path planning:

McCarthy, Zoe, Timothy Bretl, and Seth Hutchinson. "Proving path non-existence using sampling and alpha shapes."
Robotics and Automation (ICRA), 2012 IEEE International Conference on. IEEE, 2012.
