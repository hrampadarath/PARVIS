# PARVIS 
** PARVIS - Python skA Radio Visibility fIts Simulator **

A python based, command line tool to  simulate a full radio interferometric observation and export to FITSidi format, which can be imaged within CASA. This simulation was built so as to include the full Jones matrix of an
antenna i.e. the primary beam.  Currently only an offset-gregorian primary beam model is
implemented based upon the SKA1-MID dishes, although it is not difficult to include others.
See the description on primary_beam_model.txt in folder primary_beams
      
The simulator is executed as 

>python2.7 PARVIS.py  simulation.inputs
    
where doFITSsim.inputs is an input file that specifies the different 
observational parameters
Note before running this script, a new .xml file is generated via
>python2.7 createxml.py simulation.inputs

This is executed in the main python script
this is to create the FITS header file the full simulator should be executed using 

>./SimulatorMueller.sh 

here the terminal outputs are saved to a logfile in the folder "logs"
