#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, sys, time
import pandas as pd
import numpy as np

# disp, time and process warning
def disp(text):
    print('[{}] {}'.format(time.asctime(), text), file=sys.stderr)


# Mark_Reference, mark the different C-context
def Mark_Reference(ref,refmark,CG, CHG, CHH):
    print('Marking Reference Genome',file=sys.stderr)
    for cr in ref:
        # disp(cr)
        refcr, refmarkcr = ref[cr], refmark[cr]
        index = refcr.seq.find('C', 0, len(refcr)-2)
        while index >= 0:
            if refcr[index+1] == 'G': refmarkcr[index] = CG
            elif refcr[index+2] == 'G': refmarkcr[index] = CHG
            else: refmarkcr[index] = CHH
            index = refcr.seq.find('C', index+1, len(refcr)-2)
        index = refcr.seq.find('G', 2, len(refcr))
        while index >= 0:
            if refcr[index-1] == 'C': refmarkcr[index] = CG
            elif refcr[index-2] == 'C': refmarkcr[index] = CHG
            else: refmarkcr[index] = CHH
            index = refcr.seq.find('G', index+1, len(refcr))
    return refmark

