#!/usr/bin/env python

import os
from zopy.doitutils import dodict

version = "v4"

root = "python d:synth_mgx_mtx.py d:hmp1-ii_metaphlan2-mtd-qcd.pcl Stool"

experiments = [
    # no spikes
    ["null"],
    # spike seq depth only (max strength default)
    ["null-dep"           , "--spike-dep"],
    # spike bugs (max strength default)
    ["null-bug"           , "--spike-bug 0.1"],
    # spike encoding (max strength default)
    ["null-enc"           , "--spike-enc 0.1"],
    # spike expression (max strength default)
    ["true-exp"           , "--spike-exp 0.1"],
    # spike expression (med strength)
    ["true-exp-med"       , "--spike-exp 0.1", "--spike-exp-strength 0.75"],
    # spike expression (low strength)
    ["true-exp-low"       , "--spike-exp 0.1", "--spike-exp-strength 0.50"],
    # spike depth and expression simultaneously (max strength defaults)
    ["true-combo-dep-exp" , "--spike-dep --spike-exp 0.1"],
    # spike bugs and expression simultaneously (max strength defaults)
    ["true-combo-bug-exp" , "--spike-bug 0.5", "--spike-exp 0.1"],
    # spike gene group encoding (max strength default, bugs "unknown")
    ["group-null-enc"     , "--spike-enc 0.1", "--spike-groups"],
    # spike gene group expression (max strength default, bugs "unknown")
    ["group-true-exp"     , "--spike-exp 0.1", "--spike-groups"],
]

def task_run_experiment( ):
    for commands in experiments:
        basename = commands[0] = os.path.join( version, commands[0] )
        commands = [root] + commands + ["# t:{}.bug_abunds.tsv".format( basename )]
        yield dodict( commands, name=basename )
