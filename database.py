""" Database management module

 author: Matthias Wiesenberger

 The idea is that the input json file can be (uniquely) hashed to a file-name
 We thus have a hash table of simulation results The user is not supposed to
 know that we work with a hash table internally or how we store the data (other
 than we return netcdf files) and be able to select one or more entries based
 on the json input file

 Open Problems:
 - how to manage restarted simulations i.e. sims, where one input file
   generates a continuous series of netcdf files
 - how to manage simulations with more than one input file (probably the
   easiest solution is to simply concatenate the two files)

 This code may be re-inventing the wheel: NoSQL database managers may do a
 better job, e.g. the document (json) category MongoDB: manage large files
 (videos, images, PDFs, etc.) along with supporting information (metadata) that
 fits naturally into JSON documents

 On the other hand
 - we do not deal with millions of entries but a only a few hundred at best so
   a DBMS might just be overkill
 - search methods for dictionaries are offered by python as well
 - There is an additional overhead of setting up and running a mongodb server
 - the netcdf files do not really map nicely into the mongodb data model since
   files larger than 16 MB require another FileGS server to be run
 - If netcdf files are stored by the DB manager how do we get access to it
   through for example paraview
 - If json files are hidden by the database manager how do we use them to run a
   code?
"""



import json
import hashlib # for the hashing
import os.path # to check for files
import glob # for listing all files in data
import subprocess # to run the create program


def create( js):
    """Run a simulation if it does not exist yet

    For this function to work the execute.sh bash script is run
    and must contain the correct path to the Feltor library
    This function raises a subprocess.CalledProcessError error
    if the simulation crashes or the input parameters are not healthy

    Parameters:
    js (dict): the complete input file to the simulation
               see also hashinput

    Returns:
    string: filename of new entry if it did not exist before
            existing filename else

   """
    hashed = hashinput(js)
    ncfile = netcdffile( hashed)
    exists = os.path.isfile( ncfile)
    if exists:
        return ncfile
    else :
        print( "Running simulation ... ")
        #First write the json file into the database
        # so that the program can read it as input
        with open( jsonfile(hashed), 'w') as f:
            inputstring = json.dumps( js, sort_keys=True, ensure_ascii=True)
            f.write( inputstring)
            f.close()
        #Run the code to create netcdf file
        try :
            subprocess.run( ["bash", "execute.sh", jsonfile(hashed), ncfile],
                    check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            #clean up entry and escalate exception
            subprocess.run( ["rm", ncfile, jsonfile(hashed)])
            raise e
        print( " ... Done")

        return ncfile

def select( js) :
    """ Select a netcdf file based on its input parameters

    This functiont raises a ValueError exception if the file does not exist
    Parameters:
    js (dict) : The complete input file to the simulation
                see also hashinput

    Returns:
    string: filename of existing file
    """
    hashed = hashinput(js)
    ncfile = netcdffile( hashed)
    exists = os.path.isfile( ncfile)
    if not exists:
        raise ValueError( 'Entry does not exist')
    else :
        return ncfile

def table():
    """ Tabulate all exisiting data in database

    Use json.dumps(table(), indent=4) to pretty print
    Note that this table is searchable/ iteratable with standard python methods

    Returns:
    dict: {id : inputfile}
    """
    table = {}
    for filename in glob.glob('./data/*.json') :
        with open( filename, 'r') as f:
            params = json.load( f)
            table[os.path.splitext( os.path.split(filename)[1])[0]] = params
            # the key is the hash value
            f.close()
    return table

def hashinput( js):
    """Hash the input file

    Params:
    js (dict): the complete input file.

    Warning:  in order to generate a unique identifier
    js needs to be normalized in the sense that the datatype must match
    the required datatype documented (e.g. 10 in a field requiring float
    is interpreted as an integer and thus produces a different hash)

    Returns:
    string: The hexadecimal sha1 hashid of the input dictionary
    """
    inputstring = json.dumps( js, sort_keys=True, ensure_ascii=True)
    hashed = hashlib.sha1( inputstring.encode( 'utf-8') ).hexdigest()
    return hashed

def jsonfile( hashid) :
    """ Create the file path to json file from the hash id """
    return os.path.join('./data', hashid+'.json')

def netcdffile( hashid) :
    """ Create the file path to netcdf file from the hash id """
    return os.path.join('./data', hashid+'.nc')

def delete( js) :
    """ Delete an entry if it exists """
    hashed = hashinput(js)
    ncfile = netcdffile( hashed)
    exists = os.path.isfile( ncfile)
    if exists :
        subprocess.run( ["rm", ncfile, jsonfile(hashed)])

def replace( js) :
    """ Force a re-simulation: delete(js) followed by create(js) """
    delete(js)
    return create(js)


def delete_all () :
    """ Delete all data in the database """
    tab = table()
    for key in tab :
        print( key)
        delete( tab[key])
