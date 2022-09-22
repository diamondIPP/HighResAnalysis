# useful aliases for run and start the analysis
# add "source .../softdir/.bash_aliases to .bashrc
HRDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $HRDIR/venv/bin/activate
alias analyse='ipython -i $HRDIR/analyse.py -- $@'
alias converter='ipython -i $HRDIR/src/converter.py -- $@'
