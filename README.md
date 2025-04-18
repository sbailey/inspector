# DESI Data Inspector

Search targets by RA,dec,radius or by TARGETID.
Interactively view the spectra or download a FITS file with the underlying data.
The Inspector can be used interactively through a website or programatically
through URLs.

This is an exploratory work in progress and not yet ready for general use.

Stephen Bailey
Berkeley Lab
Spring 2025

## Building and running podman/docker container

```
podman build -t inspector-flask-app .
podman run -d -p 5001:5001 --name inspector inspector-flask-app
podman logs inspector
podman container prune
```


