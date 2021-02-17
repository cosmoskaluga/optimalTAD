FROM python:3.9-slim

MAINTAINER Dmitrii Smirnov "dmitrii.smirnov@phystech.edu"

COPY ./setup.py /setup.py
COPY ./config.ini /config.ini	 

WORKDIR /

RUN git clone https://github.com/Homebrew/brew ~/.linuxbrew/Homebrew \
&& mkdir ~/.linuxbrew/bin \
&& ln -s ../Homebrew/bin/brew ~/.linuxbrew/bin \
&& eval $(~/.linuxbrew/bin/brew shellenv) \
&& brew --version

RUN brew install mpich2 
RUN	brew install boost 
RUN pip3 install . 


COPY . /

ENTRYPOINT ["optimalTAD"]
