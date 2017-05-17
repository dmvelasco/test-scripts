#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-snapp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-snapp.txt
#SBATCH -J snapp
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=9000M
#SBATCH -t 10-00:00:00
set -e
set -u

# Load modules
module load java

# Declare directories
dir1="/home/dmvelasc/Software/beast/bin"			# Beast binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory

# step 1
# convert nexus format to xml with beauti - convert locally
# could probably run beauti using srun and X11, below is command and error message

#dmvelasc@farm:~/Software/beast/bin$ ./beauti
#Exception in thread "main" java.awt.HeadlessException:
#No X11 DISPLAY variable was set, but this program performed an operation which requires it.
#	at java.awt.GraphicsEnvironment.checkHeadless(GraphicsEnvironment.java:207)
#	at java.awt.Window.<init>(Window.java:536)
#	at java.awt.Frame.<init>(Frame.java:420)
#	at java.awt.Frame.<init>(Frame.java:385)
#	at beast.app.util.Utils6.startSplashScreen(Unknown Source)
#	at beast.app.beauti.BeautiLauncher.main(Unknown Source)

# step 2
# run with beast/snapp
# attempting ingroup snippet
#"$dir1"/beast "$dir2"/snapptest_outgroup_2016-10-18.xml

# works without ingroup snippet
#"$dir1"/beast "$dir2"/snapptest_2016-10-18.xml # still works, even though snapptest_amyg_2016-10-24.xml did not
#"$dir1"/beast "$dir2"/snapptest_amyg_2016-10-24.xml
# had problems with amyg_2016-10-24 due to some issues with farm
# (patching and memory allocation)
# was producing errors that seemed java related, then apache related
# comparing xml files between snapptest_2016-10-18.xml and snapptest_amyg_2016-10-24.xml did not reveal obvious differences
#"$dir1"/beast "$dir2"/amygtest_2016-10-26.xml
# works; uses the same nexus file as snapptest_amyg_2016-10-24.xml just saved with a shorter name
# the name is somewhat perpetuated in some of the java tags, so perhaps that was the issue
# below 2016-11-03
#"$dir1"/beast "$dir2"/amyg3.xml
# below 2017-03-03
"$dir1"/beast "$dir2"/snapp_amygtest4_final.xml
