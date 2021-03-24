#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse
from random import choice, random, shuffle
from collections import namedtuple

from numpy import mean, std, log
from numpy.random import normal, multinomial

import zopy.utils as zu
from zopy.table2 import table, nesteddict2table

#-------------------------------------------------------------------------------
# cli
#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser( )

parser.add_argument( "hmp" )
parser.add_argument( "site" )
parser.add_argument( "basename" )

parser.add_argument( "--n-bugs"            , type=int,   default=100 )
parser.add_argument( "--n-samples"         , type=int,   default=100 )
parser.add_argument( "--n-groups"          , type=int,   default=2000 )
parser.add_argument( "--n-pangenes"        , type=int,   default=1000 )
parser.add_argument( "--coreness"          , type=float, default=0.8 )

parser.add_argument( "--dna-read-depth"    , type=float, nargs=2, default=[10e6, 20e6] )
parser.add_argument( "--rna-read-depth"    , type=float, nargs=2, default=[ 5e6, 10e6] )
parser.add_argument( "--spike-dep"         , action="store_true" )
parser.add_argument( "--spike-dep-strength", type=float, default=1.0 )

parser.add_argument( "--spike-groups"      , action="store_true" )

parser.add_argument( "--spike-bug"         , type=float, default=0 )
parser.add_argument( "--spike-bug-strength", type=float, default=1.0 )

parser.add_argument( "--spike-enc"         , type=float, default=0 )
parser.add_argument( "--spike-enc-strength", type=float, default=1.0 )

parser.add_argument( "--spike-exp"         , type=float, default=0 )
parser.add_argument( "--spike-exp-strength", type=float, default=1.0 )

args = parser.parse_args( )

#-------------------------------------------------------------------------------
# constants
#-------------------------------------------------------------------------------

# truncation for normals in log normals
c_zmax = 3

#-------------------------------------------------------------------------------
# gene named tuple object
#-------------------------------------------------------------------------------

Gene = namedtuple( "Gene", ["bug", "group"] )

#-------------------------------------------------------------------------------
# utils
#-------------------------------------------------------------------------------

def closest_non_none( a, i ):
    # scan forward
    fore = 0
    for fore in range( 0, len( a ) - i ):
        if a[i + fore] is not None:
            break
    else:
        fore = len( a )
    # scan backward
    back = 0
    for back in range( 0, i + 1 ):
        if a[i - back] is not None:
            break
    else:
        back = len( a )
    # pick
    return (i + fore) if fore <= back else (i - back)

def spike( F, M, sign, strength, nonzero=None ):
    zeroes = set( )
    radius = 1 - strength
    if nonzero is not None:
        for d in nonzero:
            for k, v in d.items( ):
                if v == 0:
                    zeroes.add( k )
    fvals = [v for k, v in F.items( ) if k not in zeroes]
    fvals.sort( )
    mkeys = [k for k, v in M.items( ) if k not in zeroes]
    mkeys = sorted( mkeys, key=lambda x: sign * M[x] )
    for i, k in enumerate( mkeys ):
        # convert i to relative rank, r
        r = i / float( len( mkeys ) )
        # determine window bottom
        a = r - radius * r
        # determine window top
        b = r + radius * (1 - r)
        # select random position in window
        s = a + (b - a) * random( )
        # convert to integer position in fvals
        x = int( s * len( fvals ) )
        # convert to closest non-none index
        x = closest_non_none( fvals, x )
        # update spike
        F[k] = fvals[x]
        # forget selected value
        fvals[x] = None
    return None

def rename( thing ):
    ret = thing
    if type( thing ) is Gene:
        ret = "{}_{}".format( thing.bug, thing.group )
    return ret    

def write_abunds( values, depths, metadata, path, groups=False ):
    T = nesteddict2table( values )
    T.apply_rowheads( rename )
    if groups:
        T.groupby( lambda x: x.split( "_" )[1], sum )
    T.rowsort( )
    T.data[0][0] = "#"
    T.normalize_columns( )
    # colname -> sample
    for colname, col in T.iter_cols( ):
        # fixes: multinomial( ) sensitive to sum( col ) rounding to 1 + eps
        col[-1] = 1 - sum( col[:-1] )
        counts = multinomial( depths[colname], col )
        for i, c in enumerate( counts ):
            T.set( i+1, colname, int( c ) )
    M = nesteddict2table( {"Phenotype":metadata, "SeqDepth":depths} )
    T.metamerge( M )
    T.dump( path )
    return None
    
def write_spiked( d, path ):
    if len( d ) > 0:
        with open( path, "w" ) as fh:
            for item in sorted( d ):
                zu.tprint( rename( item ), d[item], file=fh )
    return None

def trunc_normal( m, sd, zmax ):
    outlier = True
    while outlier:
        sim = normal( m, sd )
        if abs( m - sim ) / sd < zmax:
            outlier = False
    return sim

#-------------------------------------------------------------------------------
# munge hmp data
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "parsing HMP data" )

T = table( args.hmp )
T.select( "STSite", args.site, transposed=True )
T.select( "VISNO", "1", transposed=True )
T.head( "SRS", invert=True )
T.apply_rowheads( lambda x: x.split( "|" )[-1] )
T.grep( "headers", "s__" )
T.grep( "headers", "_unclassified", invert=True )
T.dump( "subset.tmp" )
T.float( )
T.unrarify( 1e-20, 1 )

bugs = []
for bug, row in T.iter_rows( ):
    stats = []
    nonzero = [log( k ) / log( 10 ) for k in row if k > 0]
    stats.append( len( nonzero ) / float( len( row ) ) )
    stats.append( mean( nonzero ) )
    stats.append( std( nonzero ) )
    bugs.append( stats )

# subset top bugs
bugs = sorted( bugs, key=lambda stats: -stats[0] )[0:min( len( bugs ), args.n_bugs )]

#-------------------------------------------------------------------------------
# simulate sample names and metadata
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "simulating metadata" )

META = {}
for sample in range( args.n_samples ):
    sample = "SAMPLE{:04d}".format( sample + 1 )
    META[sample] = 0 if random( ) < 0.5 else 1

#-------------------------------------------------------------------------------
# simulate sample read depths / spike?
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "simulating sequencing depths" )

DNA_DEPTH = {}
RNA_DEPTH = {}
for deps, lims in zip( [DNA_DEPTH, RNA_DEPTH], [args.dna_read_depth, args.rna_read_depth] ):
    lo, hi = lims
    for sample in META:
        deps[sample] = int( lo + random( ) * (hi - lo) )  
    if args.spike_dep:
        spike( deps, META, 1, args.spike_dep_strength )

#-------------------------------------------------------------------------------
# simulate bugs / pangenomes
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "simulating BUG/PAN data" )

groups = ["GROUP{:06d}".format( i + 1 ) for i in range( args.n_groups )]

PAN = {}
BUG = {}
for i, stats in enumerate( bugs ):
    # simulate abundances for this bug
    bug = "BUG{:04d}".format( i + 1 )
    p, m, sd = stats
    inner = BUG.setdefault( bug, {} )
    for sample in META:
        abund = 0
        if random( ) < p:
            abund = 10 ** trunc_normal( m, sd, c_zmax )
        inner[sample] = abund
    # sample a pangenome for this bug
    shuffle( groups )
    PAN[bug] = [Gene( bug, k ) for k in sorted( groups[0:args.n_pangenes] )]
    
#-------------------------------------------------------------------------------
# spike bugs?
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "spiking bugs" )

spiked = {}
for bug in BUG:
    if random( ) < args.spike_bug:
        spiked[bug] = -1 if random( ) < 0.5 else 1
for bug, sign in spiked.items( ):
    spike( BUG[bug], META, sign, args.spike_bug_strength )

write_spiked( spiked, args.basename + ".bug_spiked.tsv" )

#-------------------------------------------------------------------------------
# compute/write bug abundance
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "writing bugs" )

write_abunds( BUG, DNA_DEPTH, META, args.basename + ".bug_abunds.tsv" )
    
#-------------------------------------------------------------------------------
# simulate gene encoding
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "simulating encoding" )

ENC = {}
for bug, bdict in BUG.items( ):
    for gene in PAN[bug]:
        inner = ENC.setdefault( gene, {} )
        for sample in bdict:
            inner[sample] = 1 if random( ) < args.coreness else 0
            
#-------------------------------------------------------------------------------
# spike gene encoding?
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "spiking encoding" )

groups = {}
if args.spike_groups:
    groups = {gene.group for gene in ENC}
    groups = {k for k in groups if random( ) < args.spike_enc}
    groups = {k:(-1 if random( ) < 0.5 else 1) for k in groups}
    
spiked = {}
for gene in ENC:
    if args.spike_groups:
        if gene.group in groups:
            spiked[gene] = groups[gene.group]
    elif random( ) < args.spike_enc:
        spiked[gene] = -1 if random( ) < 0.5 else 1
for gene, sign in spiked.items( ):
    spike( ENC[gene], META, sign, args.spike_enc_strength, nonzero=[BUG[gene.bug]] )
    
write_spiked( spiked, args.basename + ".mgx_spiked.tsv" )
write_spiked( groups, args.basename + ".mgx_spiked_groups.tsv" )

#-------------------------------------------------------------------------------
# compute/write mgx
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "computing mgx" )

MGX = {}
for gene, gdict in ENC.items( ):
    inner = MGX.setdefault( gene, {} )
    for sample in gdict:
        inner[sample] = BUG[gene.bug][sample] * gdict[sample]
    
zu.say( args.basename, "->", "writing mgx" )

write_abunds( MGX, DNA_DEPTH, META, args.basename + ".mgx_abunds.tsv" )
write_abunds( MGX, DNA_DEPTH, META, args.basename + ".mgx_abunds_groups.tsv", groups=True )
    
#-------------------------------------------------------------------------------
# simulate expression
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "simulating expression" )

# force groups to be similarly expressed over bugs
scales = {gene.group for gene in MGX}
scales = {k:trunc_normal( 0, 1, c_zmax ) for k in scales}

EXP = {}
for gene, gdict in MGX.items( ):
    inner = EXP.setdefault( gene, {} )
    scale = scales[gene.group]
    for sample in gdict:
        inner[sample] = 2 ** trunc_normal( scale, 1, c_zmax )

#-------------------------------------------------------------------------------
# spike expression?
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "spiking expression" )

groups = {}
if args.spike_groups:
    groups = {gene.group for gene in EXP}
    groups = {k for k in groups if random( ) < args.spike_exp}
    groups = {k:(-1 if random( ) < 0.5 else 1) for k in groups}

spiked = {}
for gene in EXP:
    if args.spike_groups:
        if gene.group in groups:
            spiked[gene] = groups[gene.group]
    elif random( ) < args.spike_exp:
        spiked[gene] = -1 if random( ) < 0.5 else 1
for gene, sign in spiked.items( ):
    spike( EXP[gene], META, sign, args.spike_exp_strength, nonzero=[BUG[gene.bug], ENC[gene]] )
    
write_spiked( spiked, args.basename + ".mtx_spiked.tsv" )
write_spiked( groups, args.basename + ".mtx_spiked_groups.tsv" )

#-------------------------------------------------------------------------------
# compute/write mtx
#-------------------------------------------------------------------------------

zu.say( args.basename, "->", "computing mtx" )

MTX = {}
for gene in EXP:
    inner = MTX.setdefault( gene, {} )
    for sample in EXP[gene]:
        inner[sample] = EXP[gene][sample] * MGX[gene][sample]

zu.say( args.basename, "->", "writing mtx" )

write_abunds( MTX, RNA_DEPTH, META, args.basename + ".mtx_abunds.tsv" )
write_abunds( MTX, RNA_DEPTH, META, args.basename + ".mtx_abunds_groups.tsv", groups=True )
