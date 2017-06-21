# PARVIS 

<p><b>PARVIS - Python skA Radio Visibility fIts Simulator </b></p>

A python2.7.x based, command line tool to  simulate a full radio interferometric observation and export to FITSidi format, which can be imaged within CASA. It is based upon Danny Price's pyFitsidi (https://github.com/telegraphic/pyfitsidi). An updated version is included in the "Imports" folder.

This simulation was created to investigate the effects of a time and frequency variable Jones matrix of antenna i.e. the primary beam, on the output images.  Currently only an offset-gregorian primary beam model is
implemented based upon the SKA1-MID dishes, although it is not difficult to include others.
The primary beams can be downloaded from <a href="https://www.dropbox.com/s/pmkck48tc78esda/SKA_MID_Primary_Beams.tgz?dl=0">here</a>. Please see the description on primary_beam_model.txt in folder "primary_beams"

>python2 PARVIS.py -f simulation.inputs
    
where simulation.inputs is an input file that specifies the different observational parameters
Note before running this script, a new .xml file is generated to create the FITS header file, and is executed in the main python script. 
The full simulator should be executed using 

>./SimulatorMueller.sh 

here the terminal outputs are saved to a logfile in the folder "logs"

<p><b>Inputs:</b></p>

<p><i>Models</i></p> Currently the simulator can only use point or unresolved sources. Future versions will hopefully be able to incorporate, extednded/more interesting sources. To create a list of point sources, use the point_src_list_generator.py in the "tools folder". The simulator uses an input file and is executed as follows:

<p><i>Antenna configurations</i></p> Only the VLA in C configuration and a version if the SKA1-MID baselines are included in the "antconfig" folder. However, it is possible to use any antenna configuration file. 


<p><b>Dependencies:</b></p>

numpy, scipy astropy, lxml, xml.etree.ElementTree and ephem

<p><b>Other radio interferometric astronomical simulation tools</b></p>

<ul>
  <li><a href = "https://github.com/crpurcell/friendlyVRI"> The Friendly Virtual Radio Interferometer (VRI)</a> - is a full graphic simulator written by Cormac Purcell and and Roy Truelove, from the University of Sydney.</li>
  <li><a href="https://launchpad.net/apsynsim">APSYNSIM</a>  is a 'full-fat' simulator by Ivan Marti-Vidal at Onsala Space Observatory, Sweden.</li>
  <li>

<a href="http://www.jb.man.ac.uk/pynterferometer/">Pynterferometer</a> is a public demonstration tool by Adam Avison and Sam George and has excellent documentation in the accompanying <a href="http://iopscience.iop.org/article/10.1088/0143-0807/34/1/7/meta;jsessionid=8C4EF5D281393737D6473793D4E70962.c3.iopscience.cld.iop.org">paper</a>.</li>
</ul> 



