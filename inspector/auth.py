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
    return (username == os.environ['DESI_COLLAB_USERNAME']) and (password == os.environ['DESI_COLLAB_PASSWORD'])

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

