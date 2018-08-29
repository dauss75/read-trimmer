FROM python:stretch
## read-trimmer
# Paired end trimmer for QIASeq DNA/RNA reads

# Installation of requirements
RUN pip3 install edlib Cython

# Initial installation of the trimming code
RUN git clone https://github.com/qiaseq/read-trimmer && cd read-trimmer && python3 setup.py build_ext --inplace

# install the wrapper script that will run updates before the actual trimming
RUN echo '#!/bin/bash\n\
echo "# [$(date)] Running code updates ..."\n\
cd /read-trimmer\n\
git pull --depth 1 origin master\n\
python3 setup.py build_ext --inplace\n\
echo "# [$(date)] Code updated."\n\
python3 /read-trimmer/trimmer/run.py "$@"' > /trim

RUN chmod +x /trim

ENTRYPOINT [ "/trim" ]
CMD [ "--help" ]

