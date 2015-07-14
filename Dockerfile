FROM cami/binning
# ubuntu packages
RUN apt-get update && apt-get install -y xz-utils wget python2.7 python-biopython sqlite3 libglib2.0-bin  

# PPS+
# install the tools under the /opt
ADD sw/tools /opt/tools
ADD default /dckr/etc/tasks.d/
ADD mg_all_analysis /dckr/etc/tasks.d/ 
ADD ppsp_contigs /dckr/etc/tasks.d/     
ADD s16_analysis /dckr/etc/tasks.d/ 
ADD train /dckr/etc/tasks.d/
ADD predict_contigs /dckr/etc/tasks.d/
ADD train_accuracy /dckr/etc/tasks.d/


RUN chmod u+x -R /opt/tools
# set the python variable
ENV PYTHONPATH /opt/tools/ppsplus

