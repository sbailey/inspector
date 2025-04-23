"""
Basic DESI public/private authentication

Example usage:

@app.route("/<string:specprod>/<string:foo>")
@conditional_auth
def healpix_radec_targets(specprod, foo):
    return f"Access to {specprod} with {foo=}"
"""

import os, sys
from flask import Flask, request, Response
from functools import wraps

#- Get user/pass from environment upon import to crash immediately if not set
try:
    _desi_username = os.environ['DESI_COLLAB_USERNAME']
    _desi_password = os.environ['DESI_COLLAB_PASSWORD']
except KeyError:
    print('ERROR: please set $DESI_COLLAB_USERNAME and $DESI_COLLAB_PASSWORD')
    sys.exit(1)

# Define a decorator for requiring HTTP Basic Auth
def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)
    return decorated

# Function to check if username and password are correct
def check_auth(username, password):
    global _desi_username, _desi_password
    return (username == _desi_username) and (password == _desi_password)

# Function to send a 401 response that enables basic auth
def authenticate():
    return Response(
    'Could not verify your access level for that URL.\n'
    'You have to login with proper credentials', 401,
    {'WWW-Authenticate': 'Basic realm="Login Required"'})

# Custom route decorator to conditionally apply requires_auth
def conditional_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        specprod = kwargs.get('specprod')
        if specprod not in ('fuji', 'guadalupe', 'iron', 'edr', 'dr1'):
            return requires_auth(f)(*args, **kwargs)
        return f(*args, **kwargs)
    return decorated

