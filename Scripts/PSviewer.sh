#!/bin/bash

watch -n 1 'ps -e -o pid,uname,cmd,pmem,pcpu --sort=-pcpu | head -15'
