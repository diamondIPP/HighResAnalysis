# useful aliases for run and start the analysis
# add "source .../softdir/.bash_aliases to .bashrc
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
alias analyse='ipython -i $DIR/src/dut_analysis.py -- $@'
alias make_run_plan='$DIR/src/spreadsheet.py -- $@'
# source /home/micha/software/root6/build/bin/thisroot.sh
# source /home/micha/python/py3/bin/activate
