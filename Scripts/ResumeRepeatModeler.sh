#!/bin/bash

RepeatModeler -database $1 -recoverDir $2 >& run.out
/home/drew/scripts/SendText "Finished recovering $1"
