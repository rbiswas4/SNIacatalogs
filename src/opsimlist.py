#!/usr/bin/env python

import lsst.sims.maf.db as db

garage = os.getenv('OpSimData')
print garage
dbAddress = 'sqlite:///' + garage + '/opsimblitz2_1018_sqlite.db'
oo = db.OpsimDatabase(dbAddress)
