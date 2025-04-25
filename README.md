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
podman build -t desi-inspector .
podman rm inspector
podman run -d -p 5001:5001 -v $DESI_ROOT:/desi:ro --name inspector desi-inspector
podman run -d -p 5001:5001 -v $DESI_ROOT:/desi:ro --name inspector -e DESI_COLLAB_USERNAME=... -e DESI_COLLAB_PASSWORD=... desi-inspector
podman logs inspector

http://0.0.0.0:5001
```
