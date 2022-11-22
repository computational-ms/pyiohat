pyProtista
==========

|build-status-gh| |pypi|

.. |pypi| image:: https://badge.fury.io/py/pyProtista.svg
    :target: https://badge.fury.io/py/pyProtista

.. |build-status-gh| image:: https://github.com/computational-ms/pyProtista/actions/workflows/tox_ci.yml/badge.svg
    :target: https://github.com/computational-ms/pyProtista/actions
    
.. |versions| image:: https://img.shields.io/pypi/pyversions/pyProtista.svg
        :target: https://pypi.python.org/pypi/pyProtista/


Buttom up proteomics tools come with many different output formats and while some standards, like mzIdentML, mzQuantML or mzTab exist, 
those standards have not be consistently implemented by all tools. As a result we have developed a python 
library that can be used to harmonize the outputs from identification and quantification tools for buttom up proteomics.

Currently supported 
 - Database search engines
 
   - comet 2020 01 04
   - mascot 2.6.2
   - msamanda 2
   - msfragger 3
   - msgf+ 2021.03.22
   - omssa 2.1.9
   - x!tandem alanine

 - Quantification engines
 
   - FlashLFQ
