# read-trimmer
Paired end trimmer for QIASeq DNA/RNA reads 

## Installation
```
pip3 install edlib
pip3 install Cython
git clone https://github.com/reineckef/read-trimmer.git
cd read-trimmer
python3 setup.py build_ext --inplace
```

## Usage Instructions

```
Run as python3 read-trimmer/trimmer/run.py

Usage instructions : python3 read-trimmer/trimmer/run.py --help
```

## Docker

To create a docker image just clone this repository (see above) 
and run:

```
docker build -t qiaseq/read-trimmer .
```

And then, run the trimming code inside the container:

```
docker run qiaseq/read-trimmer [ options ]
```
