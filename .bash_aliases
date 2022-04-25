# useful aliases for run and start the analysis
# add "source .../softdir/.bash_aliases to .bashrc
HRDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

alias analyse='ipython -i $HRDIR/analyse.py -- $@'
alias converter='ipython -i $HRDIR/src/converter.py -- $@'
alias raw='ipython -i $HRDIR/src/raw.py -- $@'

alias make_run_plan='$HRDIR/src/spreadsheet.py -- $@'
