######## ######## ######## ######## ######## ######## ######## ########
#
# This part fixes install_data to put data files in the same directory
# with the python library files, which is where our packages want
# them.
#
# This is essentially "smart_install_data" as used in the old
# setup.py files, except that it also understands wildcards
# and os-specific paths.  This means the module author can
# ask for data files with
#       "data/generic/*"
# instead of
#       glob.glob(os.path.join('data', 'generic', '*'))


import os
import glob

import distutils.util

import distutils.command.install_data

o =  distutils.command.install_data.install_data

# same trick as smart_install_data used: save the old run() method and
# insert our own run method ahead of it

o.old_run = o.run

def new_run ( self ) :
        """
        Hack for distutils to cause install_data to be in the same directory
        as the python library files.  Our packages expect this.
        """

        # We want our data files in the directory with the library files
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')


        # self.data_files is a list of
        #       ( destination_directory, [ source_file, source_file, source_file ] )
        #
        # We want to do wildcard expansion on all the file names.
        #
        l = [ ]
        for f in self.data_files :
            ( dest_dir, files ) = f
            fl = [ ]
            for ff in files :
                ff = distutils.util.convert_path(ff)
                ff = glob.glob(ff)
                fl.extend(ff)
            dest_dir = distutils.util.convert_path(dest_dir)
            l.append( ( dest_dir, fl ) )
        self.data_files = l

        # now use the original run() function to finish
        return distutils.command.install_data.install_data.old_run(self)

o.run = new_run
