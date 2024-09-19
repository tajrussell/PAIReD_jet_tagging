#!/bin/bash

condor_submit submit_XtoHH_500-1000_production.sub
condor_submit submit_XtoHH_1000-4000_production.sub
condor_submit submit_DY_production.sub
condor_submit submit_TT_production.sub
