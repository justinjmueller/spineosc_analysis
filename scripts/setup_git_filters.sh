#!/bin/sh
# Registers the notebook-output filter named in .gitattributes (local git config, not versioned).
git config filter.nbstrip.clean "python3 scripts/nbstrip_filter.py"
git config filter.nbstrip.smudge cat
git config filter.nbstrip.required false
echo "nbstrip filter registered for this clone"
