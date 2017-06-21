# PARVIS 
** PARVIS - Python skA Radio Visibility fIts Simulator 

A python2 based, command line tool to  simulate a full radio interferometric observation and export to FITSidi format, which can be imaged within CASA. It is based upon Danny Price's pyFitsidi (https://github.com/telegraphic/pyfitsidi). An updated version is included in the Imports folder.

This simulation was created to investigate the effects of a time and frequency variable Jones matrix of antenna i.e. the primary beam, on the output images.  Currently only an offset-gregorian primary beam model is
implemented based upon the SKA1-MID dishes, although it is not difficult to include others.
The primary beams can be downloaded from: https://www.dropbox.com/s/pmkck48tc78esda/SKA_MID_Primary_Beams.tgz?dl=0
See the description on primary_beam_model.txt in folder primary_beams
      
The simulator uses an input file and is is executed as 

>python2 PARVIS.py -f simulation.inputs
    
where simulation.inputs is an input file that specifies the different observational parameters
Note before running this script, a new .xml file is generated to create the FITS header file, and is executed in the main python script. 
The full simulator should be executed using 

>./SimulatorMueller.sh 

here the terminal outputs are saved to a logfile in the folder "logs"

** Dependencies:

numpy, scipy astropy, lxml, xml.etree.ElementTree and ephem

**Other radio interferometric astronomical simulation tools

<a href = "https://github.com/crpurcell/friendlyVRI"> The Friendly Virtual Radio Interferometer (VRI)</a> - is a full graphic simulator written by Cormac Purcell and and Roy Truelove, from the University of Sydney.

<a href="https://launchpad.net/apsynsim">APSYNSIM</a>  is a 'full-fat' simulator by Ivan Marti-Vidal at Onsala Space Observatory, Sweden.

<a href="http://www.jb.man.ac.uk/pynterferometer/">Pynterferometer</a> is a public demonstration tool by Adam Avison and Sam George and has excellent documentation in the accompanying <a href="http://iopscience.iop.org/article/10.1088/0143-0807/34/1/7/meta;jsessionid=8C4EF5D281393737D6473793D4E70962.c3.iopscience.cld.iop.org">paper</a>.

