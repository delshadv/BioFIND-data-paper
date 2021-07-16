Before running the scripts on DPUK portal you need to read below instructions:

Instructions for running Matlab using Linux with a graphical interface running from the VDI

Visit https://portal.dpuk.ukserp.ac.uk/ if you have already applied for data in https://portal.dementiasplatform.uk/apply and gotten your credentials.
In "VMware horizon client" open Windows 10 (third option is recommended)

1)	Type Xming in the search bar in the bottom left hand corner. Click Xming and it will open (there will be an X icon in the task bar)
2)	Open putty on the desktop or if you do not have putty on the desktop open it from r:/
3)	In the putty configuration window type “biofind” in the host name
4)	Change connection type to SSH
5)	In the menu to the left of the puTTY configurations screen Select SSH and click on X11 and then tick enable x11 forwarding.
6)	You can now select Session in the left hand menu and in “saved session name” you can type the name of your session i.e. BioFind.
7)	You can now open the session.
8)	From now on when you open puTTY you can go to session and open your saved session i.e. BioFind
9)	This will now take you to the Linux server where it will ask you for a Login which will be your portal credentials.
When logging in for the first time it will take a few minutes.
Because X11 has been ticked and Xming is running in the background if you open an application which has a user interface it will come up on the screen
So if we type matlab and press enter, then matlab Windows User interface will open (it may take a little while). The software Matlab is running on the Linux server but the keyboard, screen and mouse are running on your virtual desktop. So when you are accessing files you will be accessing files on the linux server not on your virtual desktop.
Warning – If you close the Linux server window you will close Matlab.
