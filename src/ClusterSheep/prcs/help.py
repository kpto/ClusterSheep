# -*- coding: utf-8 -*-
"""
Created on 03:48:37 02/12/2017

Author: Paul TO, Ka Po
Contact: kpto@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Division of Biomedical Engineering

Supervisor: Prof. Henry LAM, H. N.
Contact: kehlam@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Department of Chemical and Biomolecular Engineering

Desciption of this module:
"""

# ====BEGIN OF MODULE IMPORT====
from ClusterSheep.property import *
from ClusterSheep.prcs.template import create_template_file
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
MASCOT = """
                                =MM                                         
                              ZM= .+M,                                      
                            MM.      +M7                                    
                 MM.     .MM.           MM     ,MMM                         
                 .M ?M  MM.               MM 8M. M                          
                  M   MN                    MM.  M                          
                  M=M:                        8M.M.                         
                 MM                             8M=                         
              =MN                                  MM.                      
            MM.                                      MM.                    
          MM                                           MM.                  
        MM8.                                            .MM                 
           MMMMMMM,                               .  MMMMMMN               
                    ?MMM                     ..MMMMM                        
      /''\            M                         M            /''\         
     |    |          M.       MM       MM       MM          |    |     
      '..' MM.     .MM                           M       .MM '..'   
              MM.  ~M.                           7M   .MM 
                 MM+M                             M?MM                     
                   MM,         M.     .M           M                         
                   M             MMMMM             M.                         
                   M$                              M.                         
                   M.                              8M                         
                   M                               .M                         
                   M                                M                         
                   M                                M                        
                   M                                M                        
                   M                 $$             M                        
                .MMMMMM     MMMMM  MM  .M~ MMM,MMMMM$MMMM.                  
         .MMMMMM.      .MMMM    .OM,    .MM.      .M$   .MM                 
         M   .M           M       .                .        M                 
         MM                                                 .M                 
          M    .                                         ZM.                
           8MMM                                             MM                
              .M    MM                                     M.                
              M MMMM.MM     .M       M            MO    .MMM                
              M     8MMM. 8MMMM     +MM     M8   .M 7MMMM..               
            .M=    M .M...M M: :MMMM~  MMMMM.MMMMM= .M.M$                
            MN  .MM. MM   M+M    :MM   M.M   MDM  M MM. M:              
           .M .MM    M   ?M.M.   M.M: M7 M   M M. M  MM.M.             
            MM.     M.   M  MN  ,M  M M   M..M  M M    MM             
                    M.  M$  .M  .M  IMM    M8    MM           
                     MMM     M. ~M   $      M     M           
                     .7       MMM                                
                               M                                 
"""

HEADER = '''
{name} is a GPU accelerated program for MS2 spectral clustering.


Usage:

    {namel} [--switch-type-option --input-type-option=value] files


Examples:

    {namel} --file-list=/path/to/list --config=/path/to/config --name=run --ignore-errors
    {namel} --config=/path/to/config /path/to/files/
    {namel} --load-session=/path/to/a/session{sessext}
    {namel} --load-session=/path/to/a/session --print-session
    {namel} --list-gpus
    {namel} (shows this help and generates a configuration file template)


'''

OPTIONS = [
    ['--version', ('Print the version of this software.',)],
    ['--name=', ('Give the new session a name.',
                 'If not provided, the new session will be named by date and time.')],
    ['--fork=', ('Fork a session with a new name.',
                 'It is usually used to run clustering on prepared material for multiple times with different configuration.'),],
    ['--file-list=', ('Specify a file listing MS experiment files and identification files.',
                      'Lines starting with # will be skipped.')],
    ['--force', ('Basically try to have a complete run.',
                 'It ignores all non-critical errors and fallback to CPU if GPU failed.')],
    ['--use-cpu', ('Use CPU for clustering instead of GPU.',)],
    ['--ignore-errors', ('Ignores all non-critical errors.',)],
    ['--config=', ('Specify a file listing parameters and assigned values.',
                   'A template named "config_template" has already been created under the working directory.',
                   'Lines starting with # will be skipped.',
                   'Default parameters are used for a new session if no config file is specified.',
                   'If it is used with --load-session, the embedded config will be modified')],
    ['--dev-mode', ('Turns on developer mode.',)],
    ['--preparation-only', ('Stops after internal index building and rank transformation.',)],
    ['--stay-interactive', ('Go into cluster viewer after all processes.',)],
    ['--load-session=', ('Loads a session. File extension can be automatically filled.',)],
    ['--no-saving', ('Runs a full run and exports all clusters to a text file and deletes all other files.',)],
    ['--keep-trash', ('Keeps all intermediate files produced during the process.',)],
    ['--re-process', ('Abandons all files of a session.',)],
    ['--re-cluster', ('Abandons clustering results only.',)],
    ['--rebuild-iden-lut', ('Abandons the current identification lookup table and rebuild a new one.',
                            'User must specify identification files.')],
    ['--list-gpus', ('List all equipped CUDA device and ends program.',)],
    ['--print-session', ('Print the information of a session, must be used with --load-session.',)]
]
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def print_help():
    create_template_file()
    print('\n'.join(map(str.rstrip, MASCOT.split('\n'))))
    print(HEADER.format(name=NAME, namel=NAME.lower(), sessext=FILE_EXTENSION_SESSION))
    print('Options:\n')
    for o in OPTIONS:
        print('    {:25}'.format(o[0]), end='')
        for n, i in enumerate(o[1]):
            if n == 0: print(i)
            else: print(' '*29, i, sep='')
        print()
    print()
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
