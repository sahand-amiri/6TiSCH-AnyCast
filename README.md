# 6TiSCH-AnyCast
## Scope

6TiSCH is an IETF standardization working group that defines a complete protocol stack for ultra reliable ultra low-power wireless mesh networks.
This simulator implements the 6TiSCH protocol stack, exactly as it is standardized.
It allows you to measure the performance of a 6TiSCH network under different conditions.

Simulated protocol stack

|                                                                                                              |                                             |
|--------------------------------------------------------------------------------------------------------------|---------------------------------------------|
| [RFC6550](https://tools.ietf.org/html/rfc6550), [RFC6552](https://tools.ietf.org/html/rfc6552)               | RPL, non-storing mode, OF0                  |
| [RFC6206](https://tools.ietf.org/html/rfc6206)                                                               | Trickle Algorithm                           |
| [draft-ietf-6lo-minimal-fragment-07](https://tools.ietf.org/html/draft-ietf-6lo-minimal-fragment-07)         | 6LoWPAN Fragment Forwarding                 |
| [RFC6282](https://tools.ietf.org/html/rfc6282), [RFC4944](https://tools.ietf.org/html/rfc4944)               | 6LoWPAN Fragmentation                       |
| [draft-ietf-6tisch-msf-10](https://tools.ietf.org/html/draft-ietf-6tisch-msf-10)                             | 6TiSCH Minimal Scheduling Function (MSF)    |
| [draft-ietf-6tisch-minimal-security-15](https://tools.ietf.org/html/draft-ietf-6tisch-minimal-security-15)   | Constrained Join Protocol (CoJP) for 6TiSCH |
| [RFC8480](https://tools.ietf.org/html/rfc8480)                                                               | 6TiSCH 6top Protocol (6P)                   |
| [RFC8180](https://tools.ietf.org/html/rfc8180)                                                               | Minimal 6TiSCH Configuration                |
| [IEEE802.15.4-2015](https://ieeexplore.ieee.org/document/7460875/)                                           | IEEE802.15.4 TSCH                           |

* connectivity models
    * Pister-hack
    * k7: trace-based connectivity
* miscellaneous
    * Energy Consumption model taken from
        * [A Realistic Energy Consumption Model for TSCH Networks](http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6627960&url=http%3A%2F%2Fieeexplore.ieee.org%2Fiel7%2F7361%2F4427201%2F06627960.pdf%3Farnumber%3D6627960). Xavier Vilajosana, Qin Wang, Fabien Chraim, Thomas Watteyne, Tengfei Chang, Kris Pister. IEEE Sensors, Vol. 14, No. 2, February 2014.

## Installation

* Install Python 2.7 (or Python 3)
* Clone or download this repository
* To plot the graphs, you need Matplotlib and scipy. On Windows, Anaconda (http://continuum.io/downloads) is a good one-stop-shop.

While 6TiSCH Simulator has been tested with Python 2.7, it should work with Python 3 as well.

## Getting Started

1. Download the code:
   ```
   $ git clone https://bitbucket.org/6tisch/simulator.git
   ```
1. Install the Python dependencies:
   `cd simulator` and `pip install -r requirements.txt`
1. Execute `runSim.py` or start the GUI:
    * runSim.py
       ```
       $ cd bin
       $ python runSim.py
       ```
        * a new directory having the timestamp value as its name is created under
          `bin/simData/` (e.g., `bin/simData/20181203-161254-775`)
        * raw output data and raw charts are stored in the newly created directory
