The mimu_gnss1.mat file contains the raw data, and an interpolation algorithm was used to fill the gaps between GNSS and IMU.  
For more details, please refer to formulas 3.53 and 3.54 in the article "Improving Accuracy of Low-Cost INS/GPS for Land Applications."  
<img width="509" alt="1722265762546" src="https://github.com/user-attachments/assets/dbdab056-ea4a-4b2f-93ce-0b1b407ce9af">  
In addition to interpolation, GNSS values can also be assigned using the nearest IMU moments or the midpoint of IMU timestamps.  
The file integrated_navigation.m is a GNSS/INS integrated navigation framework with in-motion initial alignment.  
The proposing the alignment algorithm（OBA, initial alignment algorithm） are "Velocity/Position Integration Formula Part I: Application to In-Flight Coarse Alignment" and "Velocity/Position Integration Formula Part II: Application to Strapdown Inertial Navigation Computation."  
The file P.m is a framework for integrated navigation using only GNSS positioning (without initial alignment).  
The file Alignment.m contains an in-motion initial alignment algorithm for integrated navigation that I proposed.  
https://kns.cnki.net/kcms2/article/abstract?v=gisQO9UvOsZgkjGM73sYt9sN-WXA-rNXbfb276Tz4o9I-S4fvaevBwD9lXvBytq-DbNlOi-ognB62VWUIdGKEmpFRt5tDVuLgN36KEdfYavdXxjNUf3InCIC_NgpqFtfVnv4e6QzCnDmTW_ot7ZMylROpzjEjjUxrpi_1RFP6FpNVXtV5KbPnyQmHQg89Q-xRDIKA032g5Y=&uniplatform=NZKPT&language=CHS  
My knowledge is limited; this is for reference only.
