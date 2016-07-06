----- smbkits -----

A data analysis package for STM on SmB6.


____ Installation instructions ____
1. From this directory run:
   $ python setup.py develop
2. From any directory open a python prompt and run:
   > import smbkits

____ Usage ____
To import data you will need to change the 'dataFolder' variable in the script smbdata.py.  
Then you can import .3ds using:
   smbkits.importDOSmap('Name_of_data_set')
The 'tools' package contains useful functions: 
   e.g. smbkits.tools.azimuthalAverage(*args) to azimuthally average data.

