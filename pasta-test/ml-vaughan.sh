#! /bin/bash 
# Use as:
# > . ./ml-vaughan.sh [-v]
# The option -v selects listing currently loaded modules.

#----------------------------------------------------------
help() 
{
   # Display Help
   >&2 echo "Load modules for Vaughan."
   >&2 echo
   >&2 echo "Syntax : . ml-vaughan.sh [-h|--help|-v|--verbose][-t|--toolchain] [intel|gnu]]"
   >&2 echo "Options:"
   >&2 echo "  h|help      Print this Help."
   >&2 echo "  v|verbose   list modules loaded."
   >&2 echo "  t|toolchain toolchain needed 'intel'|'gnu'."
   >&2 echo
}

# Process cli arguments
OPTIONS=$(getopt -o "hvt": \
                --long "help,verbose,toolchain:" \
                 -n '$0' -- "$@")

if [ $? != 0 ] ; then echo "Failed to parse options." >&2 ; exit 1 ; fi

eval set -- "$OPTIONS" # Note the quotes around '$OPTIONS': they are essential!

HELP=false
VERBOSE=false
TOOLCHAIN=intel
while true; do
  case "$1" in
    -h | --help      ) HELP=true; shift ;;
    -v | --verbose   ) VERBOSE=true; shift ;;
    -t | --toolchain ) TOOLCHAIN="$2"; shift 2 ;;
    * ) break ;;
  esac
done

if [ "$HELP" = true ]; then help; fi

module --force purge

# the test works fine with these modules:
case "$TOOLCHAIN" in
    intel)  ml calcua/2024a;
            ml intel;
            ml iimkl;;
    gnu  )  ml calcua/2023a; 
            ml ScaLAPACK;;
    *    )  >&2 echo "Toolchain unknown: ${TOOLCHAIN}."
esac

if [ "$VERBOSE" = true ]; then 
    ml
fi