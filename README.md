# LAURA
LAURA (LetsAnalysis,Use andReport ofAsteroseismology)is written in Python 3 with 5 different packages to: download,reduce  and  analyze  raw  light  curves  under  the  Asterosismol-ogy theory. 

he code uses an object oriented approach and canbe  easily  use  even  for  non-experts  in  program  language.  Allthe  additional  python  packages  that  we  use  are  free  and  easyto  install  using  the  pip  package.  Our  code  is  summarized  inonly 2 files:LAURA_master.pyandLAURA_aux.py. From theLAURA_master.pyfile, the user can do all steps to download,reduce and analize the light curve


# INSTALLATION

This is the list of the different packages/software you will have to install in order to 
use the tool.
If you run into any problems, please contact Leandro de Almeida: monolipo.physics@gmail.com.
Any use or modification of this tool can be done, as long as the author is warned
These instructions were adapted from Morisset, C 2017.
Collaboration fellows: Bruno Lustosa, Guilherme Monteiro - UFRN 2017

In order to use our code, make sure you have python 3.x installed with the follow packages:
scipy, numpy, matplotlib, kplr, peakutils, astropy and lmfit. Installing the 
necessary packages in Ubuntu GNU/Linux is done using the follow:


    sudo apt-get install build-essential
    sudo apt-get install python3
    sudo apt-get install python3-pip
    sudo pip install numpy
    sudo pip install matplotlib
    sudo pip install scipy
    sudo pip install astropy
    sudo pip install kplr
    sudo pip install peakutils
    sudo pip install lmfit

Download our source code from github and add the path of where you want to install LAURA. In bash:

\begin{lstlisting}[language=Python]
export PYTHONPATH=/home/USER_NAME/LAURA-1.0.0/source:$PYTHONPATH
\end{lstlisting}

Add this command to your /.bashrc file or similar (/.cshrc, /.bashrc, /.profile). Go to your
LAURA folder and run the checkme.py file to verify if everything is in order. You should get 
a "good to go!" message.

git
===

You will need to have installed git software. You can verify if you
already have it by doing in a terminal: ::
   which git

if no link is given, you have to download. Depending on your operating system:

I. Mac OSX: install the package from https://git-scm.com/download/mac
II. Linux: depending on your distribution:

    A. ``sudo yum install git``
    B. ``sudo apt-get install git``

Once done, from a new terminal create a directory dedicated for
this tool with the name LAURA, and from it, clone this github files.
::   
   git clone https://github.com/monolipo/Light-Curve-Python-Analisys-Tool.git

This will download the code and instructions.

