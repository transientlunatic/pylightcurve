"""
pyLightcurve configuration file functionality
---------------------------------------------

This code is heavily based upon the configuration file
handling of the SunPy project.



"""
import os
import tempfile
import ConfigParser
import pylightcurve as lightcurve

__all__ = ['load_config', 'print_config']

def load_config():
    """
    Read the pyLightcurve configuration file. If one does not exists in the user's
    home directory then read in the defaults from module
    """
    config = ConfigParser.SafeConfigParser()

    # Get locations of pyLightcurve configuration files to be loaded
    config_files = _find_config_files()

    # Read in configuration files
    config.read(config_files)

    # Specify the working directory as a default so that the user's home
    # directory can be located in an OS-independent manner
    if not config.has_option('general', 'working_dir'):
        try:
            config.set('general', 'working_dir', os.path.join(_get_home(), "pylightcurve"))
        except ConfigParser.NoSectionError:
            # Create non-existent section
            config.add_section('general')
            config.set('general', 'working_dir', os.path.join(_get_home(), "pylightcurve"))

    if not config.has_option('downloads', 'download_dir'):
        try:
            config.set('downloads', 'download_dir', os.path.join(_get_home(), "pylightcurve"))
        except ConfigParser.NoSectionError:
            # Create non-existent section
            config.add_section('downloads')
            config.set('downloads', 'download_dir', os.path.join(_get_home(), "pylightcurve"))

    # Use absolute filepaths and adjust OS-dependent paths as needed
    filepaths = [
        ('downloads', 'download_dir')
    ]
    _fix_filepaths(config, filepaths)

    # check for pylightcurve working directory and create it if it doesn't exist
    if not os.path.isdir(config.get('downloads', 'download_dir')):
        os.mkdir(config.get('downloads', 'download_dir'))

    return config

def _is_writable_dir(p):
    """Checks to see if a directory is writable"""
    return os.path.isdir(p) and os.access(p, os.W_OK)

def _get_home():
    """Find user's home directory if possible.
    Otherwise raise error.

    """
    path = path=os.path.expanduser("~")

    if not os.path.isdir(path):
        for evar in ('HOME', 'USERPROFILE', 'TMP'):
            try:
                path = os.environ[evar]
                if os.path.isdir(path):
                    break
            except KeyError:
                pass
    if path:
        return path
    else:
        raise RuntimeError('please define environment variable $HOME')

def _find_config_files():
    """Finds locations of pyLightcurve configuration files"""
    config_files = []
    config_filename = 'pylightcurve'

    # find default configuration file
    module_dir = os.path.dirname(lightcurve.__file__)
    config_files.append(os.path.join(module_dir, 'data', 'pylightcurverc'))

    # if a user configuration file exists, add that to list of files to read
    # so that any values set there will overide ones specified in the default
    # config file
    config_path = _get_user_configdir()

    if os.path.exists(os.path.join(config_path, config_filename)):
        config_files.append(os.path.join(config_path, config_filename))

    return config_files

def _get_user_configdir():
    """
    Return the string representing the configuration dir.
    The default is "HOME/.pylightcurve".  You can override this with the
    PYLIGHTCURVE_CONFIGDIR environment variable
    """
    configdir = os.environ.get('PYLIGHTCURVE_CONFIGDIR')

    if configdir is not None:
        if not _is_writable_dir(configdir):
            raise RuntimeError('Could not write to PYLIGHTCURVE_CONFIGDIR="%s"' %
                               configdir)
        return configdir

    h = _get_home()
    p = os.path.join(_get_home(), '.pylightcurve')

    if os.path.exists(p):
        if not _is_writable_dir(p):
            raise RuntimeError("'%s' is not a writable dir; you must set %s/."
                               "pyLightcurve to be a writable dir.  You can also set "
                               "environment variable PYLIGHTCURVE_CONFIGDIR to any "
                               "writable directory where you want matplotlib "
                               "data stored " % (h, h))
    else:
        if not _is_writable_dir(h):
            raise RuntimeError("Failed to create %s/.pylightcurve; consider setting "
                               "PYLIGHTCURVE_CONFIGDIR to a writable directory for "
                               "pylightcurve configuration data" % h)

        os.mkdir(p)

    return p

def _fix_filepaths(config, filepaths):
    """Converts relative filepaths to absolute filepaths"""
    # Parse working_dir
    working_dir = _expand_filepath(config.get("general", "working_dir"))

    try:
        config.set('general', 'working_dir', working_dir)
    except ConfigParser.NoSectionError:
        # Create non-existent section
        config.add_section('general')

    for f in filepaths:
        val = config.get(*f)

        filepath = _expand_filepath(val, working_dir)

        # Create dir if it doesn't already exist
        if not os.path.isdir(filepath):
            os.makedirs(filepath)

        # Replace config value with full filepath
        params = f + (filepath,)
        try:
            config.set(*params)

        except ConfigParser.NoSectionError:
            # Create non-existent section
            config.add_section(params[0])

def _expand_filepath(filepath, working_dir=""):
    """Checks a filepath and expands it if necessary"""
    # Expand home directory
    if filepath[0] == "~":
        return os.path.abspath(os.path.expanduser(filepath))
    # Check for /tmp
    elif filepath == "/tmp":
        return tempfile.gettempdir()
    # Relative filepaths
    elif not filepath.startswith("/"):
        return os.path.join(working_dir, filepath)
    # Absolute filepath
    else:
        return filepath
