# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Update the package list and install non-python tools and dependencies
RUN apt-get update && \
    apt-get install -y make git subversion gcc libbz2-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy the current directory contents into the container at /app
COPY . /app

# Install python dependencies; final chmod is needed to be able to run a non-root
RUN pip install --no-cache-dir -r requirements.txt                  &&\
    pip install git+https://github.com/desihub/desiutil@3.5.0       &&\
    pip install git+https://github.com/desihub/desitarget@2.9.0     &&\
    pip install git+https://github.com/desihub/prospect@main        &&\
    pip install git+https://github.com/desihub/redrock@0.20.4       &&\
    install_redrock_templates                                       &&\
    pip install git+https://github.com/desihub/desispec@inventory   &&\
    chmod -R o+rX /app

# Make port 5001 available to the world outside this container
EXPOSE 5001

# Location of DESI_ROOT, must be mounted at runtime
ENV DESI_ROOT=/desi

# Location of pre-generated inventory files
ENV DESI_TARGET_INVENTORY_DIR=/inventory

# Run gunicorn server
CMD ["gunicorn", "-b", "0.0.0.0:5001", "-w", "5", "app:app"]

